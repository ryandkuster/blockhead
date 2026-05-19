import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import polars as pl

import blockhead.df_manip as dm
import blockhead.parentage as pm
import blockhead.utils as um
import blockhead.wrong_calls as wc

from blockhead.df_manip import HOM_REF, HOM_ALT, MIER_CORRECT, MIER_INCORRECT, MIER_MISSING


def _write_mier_summary(path, counts):
    with open(path, "w") as o:
        o.write("f1\tcorrect\tincorrect\tunknown\tcorrect_known\tcorrect_total\n")
        for f1, (correct, incorrect, unknown) in counts.items():
            known = correct + incorrect
            total = correct + incorrect + unknown
            ck = correct / known if known > 0 else 0.0
            ct = correct / total if total > 0 else 0.0
            o.write(f"{f1}\t{correct}\t{incorrect}\t{unknown}\t{ck:.3f}\t{ct:.3f}\n")


def _compute_chrom(chrom, pos_arr, gts_matrix, sample_idx_in_matrix,
                   named_f1_dt, adv_ls, f1_sibling_crosses, args):
    """
    Core per-chromosome computation given pre-loaded pos_arr and gts_matrix.
    Worker-safe: no matplotlib, no colors_dt, no truth_df.
    Returns (chrom, partial_mier, partial_hom, passing_pos_set, extra).
    extra is (hap_data, hap_data_f1_cross) (blockmode), wrong_df_c
    (wrong_calls), or None.
    """
    mier_dict = {}
    partial_mier = {f1: [0, 0, 0] for f1 in named_f1_dt}
    partial_hom  = {f1: [0, 0, 0] for f1 in named_f1_dt}

    for f1, (p1, p2) in named_f1_dt.items():
        mier = dm.compute_mier(
            gts_matrix[:, sample_idx_in_matrix[f1]],
            gts_matrix[:, sample_idx_in_matrix[p1]],
            gts_matrix[:, sample_idx_in_matrix[p2]],
        )
        mier_dict[f1] = mier

        partial_mier[f1][0] += int(np.sum(mier == MIER_CORRECT))
        partial_mier[f1][1] += int(np.sum(mier == MIER_INCORRECT))
        partial_mier[f1][2] += int(np.sum(mier == MIER_MISSING))

        p1_col = gts_matrix[:, sample_idx_in_matrix[p1]]
        p2_col = gts_matrix[:, sample_idx_in_matrix[p2]]
        hom_mask = (
            ((p1_col == HOM_REF) & (p2_col == HOM_ALT)) |
            ((p1_col == HOM_ALT) & (p2_col == HOM_REF))
        )
        hom_mier = mier[hom_mask]
        partial_hom[f1][0] += int(np.sum(hom_mier == MIER_CORRECT))
        partial_hom[f1][1] += int(np.sum(hom_mier == MIER_INCORRECT))
        partial_hom[f1][2] += int(np.sum(hom_mier == MIER_MISSING))

    if args.wrong_calls:
        wrong_df_c = wc.compute_wrong_calls_chrom(
            chrom, pos_arr, gts_matrix, sample_idx_in_matrix, named_f1_dt
        )
        return chrom, partial_mier, partial_hom, set(), wrong_df_c

    n_f1 = len(named_f1_dt)
    mier_stack = np.stack([mier_dict[f1] for f1 in named_f1_dt], axis=1)
    pct_correct     = (mier_stack == MIER_CORRECT).sum(axis=1) / n_f1
    pct_non_missing = (mier_stack != MIER_MISSING).sum(axis=1) / n_f1
    pass_mask = (pct_correct >= args.threshold) & (pct_non_missing >= args.non_missing)
    passing_pos = {int(p) for p in pos_arr[pass_mask]}

    hap_data = None
    hap_data_f1_cross = None
    if args.blockmode:
        filtered_mier_dict = {f1: mier_dict[f1][pass_mask] for f1 in named_f1_dt}
        if len(adv_ls) > 0:
            hap_data = dm.compute_hap_data(
                adv_ls, named_f1_dt, pos_arr[pass_mask], gts_matrix[pass_mask],
                sample_idx_in_matrix, filtered_mier_dict,
            )
        if len(f1_sibling_crosses) > 0:
            hap_data_f1_cross = dm.compute_hap_data_f1_cross(
                f1_sibling_crosses, pos_arr[pass_mask], gts_matrix[pass_mask],
                sample_idx_in_matrix, filtered_mier_dict,
            )

    extra = (hap_data, hap_data_f1_cross) if args.blockmode else None
    return chrom, partial_mier, partial_hom, passing_pos, extra


def _fetch_and_compute(vcf_path, chrom, sample_ls, sample_idx_in_matrix,
                       named_f1_dt, adv_ls, f1_sibling_crosses, args):
    """
    Worker for ProcessPoolExecutor: fetches chromosome via tabix then computes.
    Must be top-level for pickling.  No matplotlib state is touched here.
    """
    pos_arr, gts_matrix = dm.fetch_chrom(vcf_path, chrom, sample_ls,
                                          min_qual=args.quality)
    if pos_arr is None:
        empty_extra = (None, None) if args.blockmode else None
        return chrom, {f1: [0, 0, 0] for f1 in named_f1_dt}, \
               {f1: [0, 0, 0] for f1 in named_f1_dt}, set(), empty_extra

    return _compute_chrom(chrom, pos_arr, gts_matrix, sample_idx_in_matrix,
                          named_f1_dt, adv_ls, f1_sibling_crosses, args)


def _plot_chrom(args, adv_ls, named_f1_dt, f1_sibling_crosses,
                chrom, hap_data, hap_data_f1_cross, colors_dt, truth_df):
    """
    Plot worker: renders and saves figures for one chromosome.
    Returns only chrom (a small string) — no DataFrames cross the IPC pipe,
    which avoids the pipe-buffer deadlock that large pickled DataFrames cause.
    DataFrames are built in the main process via build_block_dfs().
    """
    if hap_data is not None:
        dm.plot_haplotypes(
            args, adv_ls, named_f1_dt, chrom, hap_data, colors_dt,
            truth_df=truth_df,
        )
    if hap_data_f1_cross is not None:
        dm.plot_haplotypes_f1_cross(
            args, f1_sibling_crosses, chrom, hap_data_f1_cross, colors_dt,
            truth_df=truth_df,
        )
    return chrom


def _accumulate(result, global_mier, global_hom, named_f1_dt,
                passing_positions, wrong_calls_dfs, args):
    chrom, p_mier, p_hom, p_pass, extra = result

    for f1 in named_f1_dt:
        for i in range(3):
            global_mier[f1][i] += p_mier[f1][i]
            global_hom[f1][i]  += p_hom[f1][i]

    for pos in p_pass:
        passing_positions.add((chrom, pos))

    if args.wrong_calls and extra is not None:
        wrong_calls_dfs.append(extra)


def main():
    args = um.parse_user_input()
    os.makedirs(args.outdir, exist_ok=True)

    parent_dt, cross_ls, named_f1_dt = pm.get_parentage(args)
    adv_ls = pm.get_advanced(named_f1_dt)
    f1_sibling_crosses = pm.get_f1_sibling_crosses(named_f1_dt)
    sample_ls = pm.get_sample_ls(named_f1_dt)

    from cyvcf2 import VCF as _VCF
    _v = _VCF(args.vcf, strict_gt=False)
    vcf_sample_order = _v.samples
    chrom_list = _v.seqnames  # non-empty only when VCF is tabix-indexed
    _v.close()

    sample_ls = sorted(set(sample_ls), key=lambda s: vcf_sample_order.index(s))
    sample_idx_in_matrix = {s: i for i, s in enumerate(sample_ls)}

    global_mier = {f1: [0, 0, 0] for f1 in named_f1_dt}
    global_hom  = {f1: [0, 0, 0] for f1 in named_f1_dt}
    passing_positions = set()
    out_block_df = pl.DataFrame()
    out_wrong_df = pl.DataFrame()
    wrong_calls_dfs = []

    colors_dt = dm.haplotype_colors(args) if args.blockmode else {}
    truth_df = None
    if args.blockmode and args.assess:
        truth_df = pl.read_csv(args.assess, has_header=True, separator="\t")

    do_block          = args.blockmode and not args.wrong_calls and len(adv_ls) > 0
    do_block_f1_cross = args.blockmode and not args.wrong_calls and len(f1_sibling_crosses) > 0
    do_any_block      = do_block or do_block_f1_cross
    n_chroms = max(len(chrom_list), 1) if chrom_list else 1
    n_plot_workers = min(args.threads, n_chroms)

    # hap_data kept in main process for DataFrame building after plot workers finish.
    hap_data_by_chrom = {}

    # -----------------------------------------------------------------------
    # Parallel path: one VCF worker per chromosome (requires tabix index).
    # Plot workers run concurrently — figures render while VCF I/O continues.
    # Only chrom (a string) is returned from plot workers; no IPC pipe bloat.
    # -----------------------------------------------------------------------
    if args.threads > 1 and chrom_list:
        n_vcf_workers = min(args.threads, len(chrom_list))
        print(f"parallel mode: {len(chrom_list)} chromosomes, "
              f"{n_vcf_workers} VCF workers, {n_plot_workers} plot workers",
              flush=True)

        with ProcessPoolExecutor(max_workers=n_vcf_workers) as vcf_pool, \
             ProcessPoolExecutor(max_workers=n_plot_workers) as plot_pool:

            vcf_futures = {
                vcf_pool.submit(
                    _fetch_and_compute,
                    args.vcf, chrom, sample_ls, sample_idx_in_matrix,
                    named_f1_dt, adv_ls, f1_sibling_crosses, args,
                ): chrom
                for chrom in chrom_list
            }
            plot_futures = {}

            for vcf_future in as_completed(vcf_futures):
                chrom = vcf_futures[vcf_future]
                try:
                    result = vcf_future.result()
                except Exception as e:
                    print(f"ERROR processing {chrom}: {e}", flush=True)
                    raise
                print(f"completed VCF {chrom}", flush=True)
                _accumulate(result, global_mier, global_hom, named_f1_dt,
                            passing_positions, wrong_calls_dfs, args)

                if do_any_block:
                    hap_data, hap_data_f1_cross = result[4]
                    if hap_data is not None or hap_data_f1_cross is not None:
                        hap_data_by_chrom[chrom] = (hap_data, hap_data_f1_cross)
                        pf = plot_pool.submit(
                            _plot_chrom,
                            args, adv_ls, named_f1_dt, f1_sibling_crosses,
                            chrom, hap_data, hap_data_f1_cross, colors_dt, truth_df,
                        )
                        plot_futures[pf] = chrom

            for pf in as_completed(plot_futures):
                chrom = plot_futures[pf]
                try:
                    pf.result()
                except Exception as e:
                    print(f"ERROR plotting {chrom}: {e}", flush=True)
                    raise
                print(f"completed plot {chrom}", flush=True)

    # -----------------------------------------------------------------------
    # Sequential path: single streaming pass (no tabix index required).
    # Plot workers still run concurrently — rendering overlaps with VCF streaming.
    # -----------------------------------------------------------------------
    else:
        if args.threads > 1:
            print("warning: VCF has no tabix index, falling back to sequential streaming",
                  flush=True)

        plot_futures = {}
        plot_pool_ctx = ProcessPoolExecutor(max_workers=n_plot_workers) \
                        if do_any_block else None

        try:
            for chrom, pos_arr, gts_matrix in dm.stream_chromosomes(
                    args.vcf, sample_ls, threads=args.threads,
                    min_qual=args.quality):
                result = _compute_chrom(
                    chrom, pos_arr, gts_matrix, sample_idx_in_matrix,
                    named_f1_dt, adv_ls, f1_sibling_crosses, args,
                )
                _accumulate(result, global_mier, global_hom, named_f1_dt,
                            passing_positions, wrong_calls_dfs, args)

                if do_any_block and plot_pool_ctx is not None:
                    hap_data, hap_data_f1_cross = result[4]
                    if hap_data is not None or hap_data_f1_cross is not None:
                        hap_data_by_chrom[chrom] = (hap_data, hap_data_f1_cross)
                        pf = plot_pool_ctx.submit(
                            _plot_chrom,
                            args, adv_ls, named_f1_dt, f1_sibling_crosses,
                            chrom, hap_data, hap_data_f1_cross, colors_dt, truth_df,
                        )
                        plot_futures[pf] = chrom

            for pf in as_completed(plot_futures):
                chrom = plot_futures[pf]
                try:
                    pf.result()
                except Exception as e:
                    print(f"ERROR plotting {chrom}: {e}", flush=True)
                    raise
                print(f"completed plot {chrom}", flush=True)

        finally:
            if plot_pool_ctx is not None:
                plot_pool_ctx.shutdown(wait=True)

    # -----------------------------------------------------------------------
    # Build smooth/assess DataFrames in main process from stored hap_data.
    # Fast numpy ops — no large data ever crosses the IPC pipe boundary.
    # -----------------------------------------------------------------------
    out_block_df_f1_cross = pl.DataFrame()
    out_wrong_df_f1_cross = pl.DataFrame()

    if (do_block or do_block_f1_cross) and (args.smooth or args.assess):
        for chrom, (hap_data, hap_data_f1_cross) in hap_data_by_chrom.items():
            if do_block and hap_data is not None:
                block_df_c, wrong_df_c = dm.build_block_dfs(
                    args, adv_ls, named_f1_dt, chrom, hap_data, truth_df=truth_df,
                )
                if block_df_c is not None:
                    out_block_df = pl.concat([out_block_df, block_df_c]) \
                                   if not out_block_df.is_empty() else block_df_c
                if wrong_df_c is not None:
                    out_wrong_df = pl.concat([out_wrong_df, wrong_df_c]) \
                                   if not out_wrong_df.is_empty() else wrong_df_c

            if do_block_f1_cross and hap_data_f1_cross is not None:
                block_df_c_fc, wrong_df_c_fc = dm.build_block_dfs_f1_cross(
                    args, f1_sibling_crosses, chrom, hap_data_f1_cross,
                    truth_df=truth_df,
                )
                if block_df_c_fc is not None:
                    out_block_df_f1_cross = pl.concat(
                        [out_block_df_f1_cross, block_df_c_fc]) \
                        if not out_block_df_f1_cross.is_empty() else block_df_c_fc
                if wrong_df_c_fc is not None:
                    out_wrong_df_f1_cross = pl.concat(
                        [out_wrong_df_f1_cross, wrong_df_c_fc]) \
                        if not out_wrong_df_f1_cross.is_empty() else wrong_df_c_fc

    # -----------------------------------------------------------------------
    # Post-processing outputs
    # -----------------------------------------------------------------------
    print("writing MIER summaries", flush=True)
    _write_mier_summary(os.path.join(args.outdir, "MIER_summary.tsv"), global_mier)
    _write_mier_summary(
        os.path.join(args.outdir, "MIER_summary_homozygous_parents.tsv"), global_hom)

    if args.wrong_calls:
        full_wc_df = pl.concat(wrong_calls_dfs).sort(["#CHROM", "POS"])
        mier_ls = [f"MIER_{f1}" for f1 in named_f1_dt]
        full_wc_df.write_csv(
            os.path.join(args.outdir, "MIER_per_snp_summary.tsv"),
            separator="\t", include_header=True
        )
        wc.plot_wrong_blocks(args, mier_ls, full_wc_df)
        return

    if args.outvcf:
        dm.write_filtered_vcf(args, passing_positions)

    if args.blockmode and args.smooth and not out_block_df.is_empty():
        out_block_df.sort(["CHROM", "POS"]).write_csv(
            os.path.join(args.outdir,
                         f"{len(adv_ls)}_haplotypes_{args.smooth}_smooth_blocks.tsv"),
            separator="\t", include_header=True
        )

    if args.blockmode and args.smooth and not out_block_df_f1_cross.is_empty():
        f1x_n = len(f1_sibling_crosses)
        out_block_df_f1_cross.sort(["CHROM", "POS"]).write_csv(
            os.path.join(args.outdir,
                         f"{f1x_n}_f1_cross_haplotypes_{args.smooth}_smooth_blocks.tsv"),
            separator="\t", include_header=True
        )

    if args.blockmode and args.assess and not out_wrong_df.is_empty():
        out_wrong_df = out_wrong_df.filter(
            pl.any_horizontal([pl.col(c).is_not_null() for c in adv_ls])
        ).sort(["CHROM", "POS"])
        out_wrong_df.write_csv(
            os.path.join(args.outdir,
                         f"{len(adv_ls)}_haplotypes_assess_blocks.tsv"),
            separator="\t", include_header=True
        )

    if args.blockmode and args.assess and not out_wrong_df_f1_cross.is_empty():
        f1x_children = [d["child"] for d in f1_sibling_crosses]
        out_wrong_df_f1_cross = out_wrong_df_f1_cross.filter(
            pl.any_horizontal([pl.col(c).is_not_null() for c in f1x_children])
        ).sort(["CHROM", "POS"])
        out_wrong_df_f1_cross.write_csv(
            os.path.join(args.outdir,
                         f"{len(f1_sibling_crosses)}_f1_cross_haplotypes_assess_blocks.tsv"),
            separator="\t", include_header=True
        )


if __name__ == "__main__":
    main()
