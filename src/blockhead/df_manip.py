import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from blockhead.parentage import get_advanced_lineage

# Genotype encoding (cyvcf2 gt_types)
HOM_REF = 0
HET = 1
UNKNOWN = 2
HOM_ALT = 3

# MIER result encoding
MIER_CORRECT = 1
MIER_INCORRECT = 2
MIER_MISSING = 3


def fetch_chrom(vcf_path, chrom, sample_names, threads=1, min_qual=None):
    """
    Fetch a single chromosome from a tabix-indexed VCF.
    Returns (pos_arr, gts_matrix) or (None, None) if no biallelic SNPs found.
    """
    from cyvcf2 import VCF
    vcf = VCF(vcf_path, strict_gt=False, threads=threads)
    vcf_samples = vcf.samples
    try:
        sample_indices = np.array([vcf_samples.index(s) for s in sample_names])
    except ValueError as e:
        raise ValueError(f"Sample not found in VCF: {e}")

    n_samples = len(sample_indices)
    _BUF = 131_072
    pos_buf = np.empty(_BUF, dtype=np.int32)
    gts_buf = np.empty((_BUF, n_samples), dtype=np.uint8)
    n = 0

    for variant in vcf(chrom):
        # if (len(variant.REF) != 1 or len(variant.ALT) != 1
        #         or len(variant.ALT[0]) != 1):
        #     continue
        if len(variant.ALT) != 1:  # exclude multiallelic sites (comma in ALT field)
            continue
        if min_qual is not None and (variant.QUAL is None or variant.QUAL < min_qual):
            continue
        if n >= len(pos_buf):
            new_size = len(pos_buf) * 2
            pos_buf = np.resize(pos_buf, new_size)
            gts_buf = np.resize(gts_buf, (new_size, n_samples))
        pos_buf[n] = variant.POS
        gts_buf[n] = variant.gt_types[sample_indices]
        n += 1

    vcf.close()
    if n == 0:
        return None, None
    return pos_buf[:n].copy(), gts_buf[:n].copy()


def stream_chromosomes(vcf_path, sample_names, threads=1, min_qual=None):
    """
    Generator yielding (chrom, pos_array, gts_matrix) for each chromosome.
    gts_matrix shape: (n_variants, n_samples), dtype uint8
    Encoding: 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT  (cyvcf2 gt_types)
    Only biallelic SNPs are yielded.
    """
    from cyvcf2 import VCF
    vcf = VCF(vcf_path, strict_gt=False, threads=threads)
    vcf_samples = vcf.samples
    try:
        sample_indices = np.array([vcf_samples.index(s) for s in sample_names])
    except ValueError as e:
        raise ValueError(f"Sample not found in VCF: {e}")

    n_samples = len(sample_indices)
    _BUF = 131_072  # initial per-chromosome row buffer

    current_chrom = None
    pos_buf = None
    gts_buf = None
    n = 0
    total_variants = 0
    biallelic_count = 0

    def _flush(chrom, pos_b, gts_b, count):
        return (chrom,
                pos_b[:count].copy(),
                gts_b[:count].copy())

    print("streaming vcf", flush=True)
    for variant in vcf:
        total_variants += 1
        # if (len(variant.REF) != 1 or len(variant.ALT) != 1
        #         or len(variant.ALT[0]) != 1):
        #     continue
        if len(variant.ALT) != 1:  # exclude multiallelic sites (comma in ALT field)
            continue
        biallelic_count += 1
        if min_qual is not None and (variant.QUAL is None or variant.QUAL < min_qual):
            continue

        chrom = variant.CHROM
        if chrom != current_chrom:
            if current_chrom is not None:
                yield _flush(current_chrom, pos_buf, gts_buf, n)
            current_chrom = chrom
            pos_buf = np.empty(_BUF, dtype=np.int32)
            gts_buf = np.empty((_BUF, n_samples), dtype=np.uint8)
            n = 0
            print(f"  chromosome: {chrom}", flush=True)

        if n >= len(pos_buf):
            new_size = len(pos_buf) * 2
            pos_buf = np.resize(pos_buf, new_size)
            gts_buf = np.resize(gts_buf, (new_size, n_samples))

        pos_buf[n] = variant.POS
        gts_buf[n] = variant.gt_types[sample_indices]
        n += 1

    if current_chrom is not None:
        yield _flush(current_chrom, pos_buf, gts_buf, n)

    vcf.close()
    print(f"{total_variants} total variants, {biallelic_count} biallelic SNPs",
          flush=True)


def compute_mier(f1_col, p1_col, p2_col):
    """
    Vectorised MIER for one trio. All inputs are uint8 numpy arrays.
    Returns int8 array: 1=CORRECT, 2=INCORRECT, 3=MISSING
    """
    missing = (f1_col == UNKNOWN) | (p1_col == UNKNOWN) | (p2_col == UNKNOWN)
    correct = (
        ((f1_col == HOM_REF) & (p1_col != HOM_ALT) & (p2_col != HOM_ALT)) |
        ((f1_col == HET) &
         ~((p1_col == HOM_REF) & (p2_col == HOM_REF)) &
         ~((p1_col == HOM_ALT) & (p2_col == HOM_ALT))) |
        ((f1_col == HOM_ALT) & (p1_col != HOM_REF) & (p2_col != HOM_REF))
    )
    return np.where(missing, MIER_MISSING,
                    np.where(correct, MIER_CORRECT, MIER_INCORRECT)).astype(np.int8)


def haplotype_colors(args):
    if args.colors:
        colors_df = pl.read_csv(args.colors, separator="\t", has_header=False,
                                new_columns=["id", "color"])
        return dict(zip(colors_df["id"], colors_df["color"]))
    return {}


def infer_haplotypes(pos_arr, gts_matrix, sample_idx_in_matrix, mier_dict, lin_dt):
    """
    Filter to sites where p1/p2 are homozygous-opposite and p3 is homozygous,
    keeping only MIER-correct calls.  Returns (pos_filtered, parent_origin) as
    numpy arrays.  All operations are vectorised numpy — no DataFrame is built.
    Integer encoding: HOM_REF=0, HET=1, UNKNOWN=2, HOM_ALT=3
    """
    p1_col  = gts_matrix[:, sample_idx_in_matrix[lin_dt["p1"]]]
    p2_col  = gts_matrix[:, sample_idx_in_matrix[lin_dt["p2"]]]
    p3_col  = gts_matrix[:, sample_idx_in_matrix[lin_dt["p3"]]]
    adv_col = gts_matrix[:, sample_idx_in_matrix[lin_dt["adv"]]]
    mier_f1  = mier_dict[lin_dt["f1"]]
    mier_adv = mier_dict[lin_dt["adv"]]

    mask = (
        (((p1_col == HOM_REF) & (p2_col == HOM_ALT)) |
         ((p1_col == HOM_ALT) & (p2_col == HOM_REF))) &
        (mier_f1  == MIER_CORRECT) &
        (mier_adv == MIER_CORRECT) &
        ((p3_col == HOM_REF) | (p3_col == HOM_ALT))
    )

    pos_f  = pos_arr[mask]
    adv_f  = adv_col[mask]
    p1_f   = p1_col[mask]
    p2_f   = p2_col[mask]
    p3_f   = p3_col[mask]

    # Infer which grandparent genotype the adv allele matches
    parent_type = np.full(len(pos_f), 99, dtype=np.int8)
    parent_type[(adv_f == HOM_REF) & (p3_f == HOM_REF)] = HOM_REF
    parent_type[(adv_f == HET)     & (p3_f == HOM_REF)] = HOM_ALT
    parent_type[(adv_f == HET)     & (p3_f == HOM_ALT)] = HOM_REF
    parent_type[(adv_f == HOM_ALT) & (p3_f == HOM_ALT)] = HOM_ALT

    parent_origin = np.where(parent_type == p1_f, lin_dt["p1"],
                    np.where(parent_type == p2_f, lin_dt["p2"], "unknown"))

    return pos_f, parent_origin


def infer_haplotypes_f1_cross(pos_arr, gts_matrix, sample_idx_in_matrix,
                              mier_dict, f1_cross_dt):
    """
    Haplotype inference for F1×F1 sibling crosses.

    Informative sites require:
      - p1/p2 homozygous-opposite (one HOM_REF, one HOM_ALT)
      - both f1a and f1b heterozygous at the site
      - MIER-correct calls for both f1a and f1b
      - child homozygous (HOM_REF or HOM_ALT) — HET is ambiguous and excluded

    Parent origin is assigned by matching the child's homozygous genotype to
    whichever of p1/p2 carries the same genotype at that site.
    """
    p1_col    = gts_matrix[:, sample_idx_in_matrix[f1_cross_dt["p1"]]]
    p2_col    = gts_matrix[:, sample_idx_in_matrix[f1_cross_dt["p2"]]]
    f1a_col   = gts_matrix[:, sample_idx_in_matrix[f1_cross_dt["f1a"]]]
    f1b_col   = gts_matrix[:, sample_idx_in_matrix[f1_cross_dt["f1b"]]]
    child_col = gts_matrix[:, sample_idx_in_matrix[f1_cross_dt["child"]]]
    mier_f1a  = mier_dict[f1_cross_dt["f1a"]]
    mier_f1b  = mier_dict[f1_cross_dt["f1b"]]

    mask = (
        (((p1_col == HOM_REF) & (p2_col == HOM_ALT)) |
         ((p1_col == HOM_ALT) & (p2_col == HOM_REF))) &
        (f1a_col  == HET) &
        (f1b_col  == HET) &
        (mier_f1a == MIER_CORRECT) &
        (mier_f1b == MIER_CORRECT) &
        ((child_col == HOM_REF) | (child_col == HOM_ALT))
    )

    pos_f   = pos_arr[mask]
    child_f = child_col[mask]
    p1_f    = p1_col[mask]
    p2_f    = p2_col[mask]

    parent_origin = np.where(child_f == p1_f, f1_cross_dt["p1"],
                    np.where(child_f == p2_f, f1_cross_dt["p2"], "unknown"))

    return pos_f, parent_origin


def compute_hap_data_f1_cross(f1_sibling_crosses, pos_arr, gts_matrix,
                               sample_idx_in_matrix, mier_dict):
    """
    Like compute_hap_data but for F1×F1 sibling crosses.
    Worker-safe: no matplotlib, no Polars.

    Returns a dict:
      {
        'max_pos':      int,
        'cross_results': [(f1_cross_dt, pos_f, parent_int8), ...]
      }
    where parent_int8 encodes 0=p1, 1=p2, 2=unknown.
    """
    cross_results = []
    for f1_cross_dt in f1_sibling_crosses:
        pos_f, parent_origin = infer_haplotypes_f1_cross(
            pos_arr, gts_matrix, sample_idx_in_matrix, mier_dict, f1_cross_dt)
        parent_int8 = np.full(len(pos_f), 2, dtype=np.int8)
        parent_int8[parent_origin == f1_cross_dt["p1"]] = 0
        parent_int8[parent_origin == f1_cross_dt["p2"]] = 1
        cross_results.append((f1_cross_dt, pos_f, parent_int8))

    return {
        'max_pos':      int(pos_arr.max()) if len(pos_arr) > 0 else 0,
        'cross_results': cross_results,
    }


def compute_hap_data(adv_ls, named_f1_dt, pos_arr, gts_matrix,
                     sample_idx_in_matrix, mier_dict):
    """
    Worker-safe: run infer_haplotypes for every adv on one chromosome.
    No matplotlib, no Polars, no large string arrays cross the IPC boundary.

    Returns a dict:
      {
        'max_pos':     int,
        'adv_results': [(adv, lin_dt, pos_f, parent_int8), ...]
      }
    where parent_int8 is an int8 array with values 0=p1, 1=p2, 2=unknown.
    """
    adv_results = []
    for adv in adv_ls:
        lin_dt = get_advanced_lineage(adv, named_f1_dt)
        pos_f, parent_origin = infer_haplotypes(
            pos_arr, gts_matrix, sample_idx_in_matrix, mier_dict, lin_dt)
        parent_int8 = np.full(len(pos_f), 2, dtype=np.int8)
        parent_int8[parent_origin == lin_dt["p1"]] = 0
        parent_int8[parent_origin == lin_dt["p2"]] = 1
        adv_results.append((adv, lin_dt, pos_f, parent_int8))

    return {
        'max_pos':     int(pos_arr.max()) if len(pos_arr) > 0 else 0,
        'adv_results': adv_results,
    }


def plot_haplotypes(args, adv_ls, named_f1_dt, chrom, hap_data,
                    colors_dt, truth_df=None):
    """
    Main-process only: matplotlib + smooth/assess.  Receives compact hap_data
    dict produced by compute_hap_data() (or built inline in sequential mode).
    Returns (block_df_chrom, wrong_df_chrom); either may be None.
    """
    alpha_val = 0.5
    recomb_ls = []

    fig, axes = plt.subplots(
        nrows=len(adv_ls), ncols=3, sharex='col',
        figsize=(32, len(adv_ls) / 3),
        gridspec_kw={'hspace': 0.3, 'width_ratios': [90, 1, 1], 'wspace': 0.03},
        squeeze=False,
    )

    max_pos = hap_data['max_pos']
    x_ticks = list(range(0, max_pos, 5_000_000))
    tick_labels = [f'{x/1e6:.0f}Mb' if x != 0 else '0' for x in x_ticks]

    block_df_chrom = None
    wrong_df_chrom = None

    from matplotlib.colors import to_rgba
    rgba_cache = {name: to_rgba(col) for name, col in colors_dt.items()}

    for idx, (adv, lin_dt, pos_f, parent_int8) in enumerate(hap_data['adv_results']):
        print(f"  processing advanced hybrid: {adv}", flush=True)

        # Reconstruct string parent_origin from compact int8
        parent_origin = np.where(parent_int8 == 0, lin_dt["p1"],
                        np.where(parent_int8 == 1, lin_dt["p2"], "unknown"))

        if args.smooth:
            s_arr     = np.where(parent_origin == lin_dt["p1"], 0, 2)
            s_new     = median_filter(s_arr, args.smooth)
            y_new_arr = np.where(s_new == 0, lin_dt["p1"], lin_dt["p2"])
            process_diffs(pos_f, parent_origin, y_new_arr)
            recomb_ls += break_points(pos_f, y_new_arr)
            label_arr, c1, c2 = y_new_arr, rgba_cache[lin_dt["p1"]], rgba_cache[lin_dt["p2"]]
            plot_x = pos_f
        else:
            label_arr, c1, c2 = parent_origin, rgba_cache[lin_dt["p1"]], rgba_cache[lin_dt["p2"]]
            plot_x = pos_f

        y_zeros = np.zeros(len(plot_x), dtype=np.float32)
        mask_p1 = (label_arr == lin_dt["p1"])
        for sel, col_val in ((mask_p1, c1), (~mask_p1, c2)):
            if sel.any():
                axes[idx, 0].scatter(plot_x[sel], y_zeros[sel], s=500, marker="|",
                                     color=col_val, alpha=alpha_val, rasterized=True)
        axes[idx, 0].set_yticks([])
        axes[idx, 0].set_ylabel(f"{adv}", labelpad=100, va="center", ha="left", rotation=0)
        axes[idx, 1].set_facecolor(colors_dt[lin_dt["f1"]])
        axes[idx, 2].set_facecolor(colors_dt[lin_dt["p3"]])

    for idx in range(len(adv_ls)):
        for col in (1, 2):
            axes[idx, col].set_yticklabels([])
            axes[idx, col].set_yticks([])
            axes[idx, col].set_xticklabels([])
            axes[idx, col].set_xticks([])

    for ax in axes:
        ax[0].autoscale(enable=True, axis='x', tight=True)

    axes[0, 0].set_title(f"{chrom}", pad=20)
    axes[0, 1].set_title("F1", pad=20)
    axes[0, 2].set_title("P3", pad=20)
    axes[-1, 0].set_xticks(x_ticks)
    axes[-1, 0].set_xticklabels(tick_labels)

    if not args.smooth:
        img_path = os.path.join(args.outdir, f"{chrom}_{len(adv_ls)}_haplotypes.png")
        print(f"  saving {img_path}", flush=True)
        plt.savefig(img_path, dpi=300, bbox_inches="tight")
    else:
        img_path = os.path.join(args.outdir, f"{chrom}_{len(adv_ls)}_haplotypes_smooth.png")
        print(f"  saving {img_path}", flush=True)
        plt.savefig(img_path, dpi=300, bbox_inches="tight")

        img_path = os.path.join(args.outdir, f"{chrom}_{len(adv_ls)}_haplotypes_breaks.png")
        fig.subplots_adjust(top=0.8)
        pos0 = axes[0, 0].get_position()
        hist_ax = fig.add_axes(
            [pos0.x0, pos0.y0 + pos0.height * 1.05, pos0.width, 0.15])
        bin_no = int(max_pos // 1e6) * 2
        hist_ax.hist(recomb_ls, bins=bin_no, range=(0, max_pos), color="grey")
        hist_ax.xaxis.set_visible(False)
        hist_ax.margins(x=0)
        hist_ax.set_title(f"{chrom}", pad=20)
        print(f"  saving {img_path}", flush=True)
        plt.savefig(img_path, dpi=300, bbox_inches="tight")

    plt.close(fig)


def build_block_dfs(args, adv_ls, named_f1_dt, chrom, hap_data, truth_df=None):
    """
    Build smooth and/or assess DataFrames in the main process from hap_data.
    Called after plot workers finish so no large data crosses the IPC pipe.
    The smooth/assess numpy ops are fast; recomputing them here is cheap.
    Returns (block_df_chrom, wrong_df_chrom); either may be None.
    """
    smooth_collected = []
    assess_collected = []

    for adv, lin_dt, pos_f, parent_int8 in hap_data['adv_results']:
        parent_origin = np.where(parent_int8 == 0, lin_dt["p1"],
                        np.where(parent_int8 == 1, lin_dt["p2"], "unknown"))

        if args.smooth:
            s_arr     = np.where(parent_origin == lin_dt["p1"], 0, 2)
            s_new     = median_filter(s_arr, args.smooth)
            y_new_arr = np.where(s_new == 0, lin_dt["p1"], lin_dt["p2"])
            smooth_collected.append((pos_f, y_new_arr, adv))

        if args.assess and truth_df is not None:
            x = pos_f.tolist()
            y = parent_origin.tolist()
            adv_truth_df = truth_df.filter(pl.col(adv).is_not_null()).sort(["CHROM", "POS"])
            truth_x = adv_truth_df.filter(pl.col("CHROM") == chrom)["POS"].to_list()
            truth_y = adv_truth_df.filter(pl.col("CHROM") == chrom)[adv].to_list()
            new_x, new_y = assess_blocks(x, y, truth_x, truth_y)
            assess_collected.append((np.asarray(new_x), np.asarray(new_y), adv))

    block_df_chrom = None
    wrong_df_chrom = None

    if smooth_collected:
        from functools import reduce
        dfs = [
            pl.DataFrame({"POS": pl.Series("POS", pos_c),
                          adv_c: pl.Series(adv_c, y_c.tolist())})
            for pos_c, y_c, adv_c in smooth_collected
        ]
        merged = reduce(
            lambda a, b: a.join(b, on="POS", how="outer", coalesce=True), dfs)
        block_df_chrom = merged.with_columns(pl.lit(chrom).alias("CHROM")) \
                               .select(["CHROM", "POS"] + [c[2] for c in smooth_collected]) \
                               .sort("POS")

    if assess_collected:
        from functools import reduce
        dfs = [
            pl.DataFrame({"POS": pl.Series("POS", x_c.astype(np.int32)),
                          adv_c: pl.Series(adv_c, y_c.tolist())})
            for x_c, y_c, adv_c in assess_collected
        ]
        merged = reduce(
            lambda a, b: a.join(b, on="POS", how="outer", coalesce=True), dfs)
        wrong_df_chrom = merged.with_columns(pl.lit(chrom).alias("CHROM")) \
                               .select(["CHROM", "POS"] + [c[2] for c in assess_collected]) \
                               .sort("POS")

    return block_df_chrom, wrong_df_chrom


def plot_haplotypes_f1_cross(args, f1_sibling_crosses, chrom,
                              hap_data_f1_cross, colors_dt, truth_df=None):
    """
    Like plot_haplotypes but for F1×F1 sibling crosses.
    Columns: haplotype track | F1a color swatch | F1b color swatch.
    """
    n = len(f1_sibling_crosses)
    alpha_val = 0.5
    recomb_ls = []

    fig, axes = plt.subplots(
        nrows=n, ncols=3, sharex='col',
        figsize=(32, n / 3),
        gridspec_kw={'hspace': 0.3, 'width_ratios': [90, 1, 1], 'wspace': 0.03},
        squeeze=False,
    )

    max_pos = hap_data_f1_cross['max_pos']
    x_ticks = list(range(0, max_pos, 5_000_000))
    tick_labels = [f'{x/1e6:.0f}Mb' if x != 0 else '0' for x in x_ticks]

    from matplotlib.colors import to_rgba
    rgba_cache = {name: to_rgba(col) for name, col in colors_dt.items()}

    for idx, (f1_cross_dt, pos_f, parent_int8) in \
            enumerate(hap_data_f1_cross['cross_results']):
        child = f1_cross_dt["child"]
        print(f"  processing f1 sibling cross: {child}", flush=True)

        parent_origin = np.where(parent_int8 == 0, f1_cross_dt["p1"],
                        np.where(parent_int8 == 1, f1_cross_dt["p2"], "unknown"))

        if args.smooth:
            s_arr     = np.where(parent_origin == f1_cross_dt["p1"], 0, 2)
            s_new     = median_filter(s_arr, args.smooth)
            y_new_arr = np.where(s_new == 0, f1_cross_dt["p1"], f1_cross_dt["p2"])
            process_diffs(pos_f, parent_origin, y_new_arr)
            recomb_ls += break_points(pos_f, y_new_arr)
            label_arr = y_new_arr
            c1 = rgba_cache[f1_cross_dt["p1"]]
            c2 = rgba_cache[f1_cross_dt["p2"]]
            plot_x = pos_f
        else:
            label_arr = parent_origin
            c1 = rgba_cache[f1_cross_dt["p1"]]
            c2 = rgba_cache[f1_cross_dt["p2"]]
            plot_x = pos_f

        y_zeros = np.zeros(len(plot_x), dtype=np.float32)
        mask_p1 = (label_arr == f1_cross_dt["p1"])
        for sel, col_val in ((mask_p1, c1), (~mask_p1, c2)):
            if sel.any():
                axes[idx, 0].scatter(plot_x[sel], y_zeros[sel], s=500, marker="|",
                                     color=col_val, alpha=alpha_val, rasterized=True)
        axes[idx, 0].set_yticks([])
        axes[idx, 0].set_ylabel(f"{child}", labelpad=100, va="center",
                                ha="left", rotation=0)
        axes[idx, 1].set_facecolor(colors_dt[f1_cross_dt["f1a"]])
        axes[idx, 2].set_facecolor(colors_dt[f1_cross_dt["f1b"]])

    for idx in range(n):
        for col in (1, 2):
            axes[idx, col].set_yticklabels([])
            axes[idx, col].set_yticks([])
            axes[idx, col].set_xticklabels([])
            axes[idx, col].set_xticks([])

    for ax in axes:
        ax[0].autoscale(enable=True, axis='x', tight=True)

    axes[0, 0].set_title(f"{chrom}", pad=20)
    axes[0, 1].set_title("F1a", pad=20)
    axes[0, 2].set_title("F1b", pad=20)
    axes[-1, 0].set_xticks(x_ticks)
    axes[-1, 0].set_xticklabels(tick_labels)

    if not args.smooth:
        img_path = os.path.join(args.outdir,
                                f"{chrom}_{n}_f1_cross_haplotypes.png")
        print(f"  saving {img_path}", flush=True)
        plt.savefig(img_path, dpi=300, bbox_inches="tight")
    else:
        img_path = os.path.join(args.outdir,
                                f"{chrom}_{n}_f1_cross_haplotypes_smooth.png")
        print(f"  saving {img_path}", flush=True)
        plt.savefig(img_path, dpi=300, bbox_inches="tight")

        img_path = os.path.join(args.outdir,
                                f"{chrom}_{n}_f1_cross_haplotypes_breaks.png")
        fig.subplots_adjust(top=0.8)
        pos0 = axes[0, 0].get_position()
        hist_ax = fig.add_axes(
            [pos0.x0, pos0.y0 + pos0.height * 1.05, pos0.width, 0.15])
        bin_no = int(max_pos // 1e6) * 2
        hist_ax.hist(recomb_ls, bins=bin_no, range=(0, max_pos), color="grey")
        hist_ax.xaxis.set_visible(False)
        hist_ax.margins(x=0)
        hist_ax.set_title(f"{chrom}", pad=20)
        print(f"  saving {img_path}", flush=True)
        plt.savefig(img_path, dpi=300, bbox_inches="tight")

    plt.close(fig)


def build_block_dfs_f1_cross(args, f1_sibling_crosses, chrom,
                              hap_data_f1_cross, truth_df=None):
    """
    Like build_block_dfs but for F1×F1 sibling crosses.
    Returns (block_df_chrom, wrong_df_chrom); either may be None.
    """
    smooth_collected = []
    assess_collected = []

    for f1_cross_dt, pos_f, parent_int8 in hap_data_f1_cross['cross_results']:
        child = f1_cross_dt["child"]
        parent_origin = np.where(parent_int8 == 0, f1_cross_dt["p1"],
                        np.where(parent_int8 == 1, f1_cross_dt["p2"], "unknown"))

        if args.smooth:
            s_arr     = np.where(parent_origin == f1_cross_dt["p1"], 0, 2)
            s_new     = median_filter(s_arr, args.smooth)
            y_new_arr = np.where(s_new == 0, f1_cross_dt["p1"], f1_cross_dt["p2"])
            smooth_collected.append((pos_f, y_new_arr, child))

        if args.assess and truth_df is not None:
            x = pos_f.tolist()
            y = parent_origin.tolist()
            child_truth_df = truth_df.filter(
                pl.col(child).is_not_null()).sort(["CHROM", "POS"])
            truth_x = child_truth_df.filter(
                pl.col("CHROM") == chrom)["POS"].to_list()
            truth_y = child_truth_df.filter(
                pl.col("CHROM") == chrom)[child].to_list()
            new_x, new_y = assess_blocks(x, y, truth_x, truth_y)
            assess_collected.append((np.asarray(new_x), np.asarray(new_y), child))

    block_df_chrom = None
    wrong_df_chrom = None

    if smooth_collected:
        from functools import reduce
        dfs = [
            pl.DataFrame({"POS": pl.Series("POS", pos_c),
                          child_c: pl.Series(child_c, y_c.tolist())})
            for pos_c, y_c, child_c in smooth_collected
        ]
        merged = reduce(
            lambda a, b: a.join(b, on="POS", how="outer", coalesce=True), dfs)
        block_df_chrom = merged.with_columns(pl.lit(chrom).alias("CHROM")) \
                               .select(["CHROM", "POS"] + [c[2] for c in smooth_collected]) \
                               .sort("POS")

    if assess_collected:
        from functools import reduce
        dfs = [
            pl.DataFrame({"POS": pl.Series("POS", x_c.astype(np.int32)),
                          child_c: pl.Series(child_c, y_c.tolist())})
            for x_c, y_c, child_c in assess_collected
        ]
        merged = reduce(
            lambda a, b: a.join(b, on="POS", how="outer", coalesce=True), dfs)
        wrong_df_chrom = merged.with_columns(pl.lit(chrom).alias("CHROM")) \
                               .select(["CHROM", "POS"] + [c[2] for c in assess_collected]) \
                               .sort("POS")

    return block_df_chrom, wrong_df_chrom



    """Second streaming pass: write variants at passing (chrom, pos) positions."""
    from cyvcf2 import VCF, Writer
    vcf_in = VCF(args.vcf, strict_gt=False)
    out_path = os.path.join(args.outdir, "MIER_filtered.vcf")
    vcf_out = Writer(out_path, vcf_in)
    written = 0
    for variant in vcf_in:
        if (variant.CHROM, variant.POS) in passing_positions:
            vcf_out.write_record(variant)
            written += 1
    vcf_out.close()
    vcf_in.close()
    print(f"wrote {written} filtered variants to {out_path}", flush=True)


# ---------------------------------------------------------------------------
# Smoothing and assessment
# ---------------------------------------------------------------------------

def median_filter(s, r):
    """
    Sliding-window median filter over a list of 0/2 values.  Because the
    only possible values are 0 and 2, the median reduces to a majority vote
    computable from a cumulative sum in O(1) per position — no per-element
    Python list slicing or statistics.median() call required.

    Window-boundary rules are identical to the original implementation:
      - very close to ends (< end_dist):  fixed small window
      - near ends (< r):                  growing/shrinking symmetric window
      - interior:                         full window of width 2r+1
    A tie (equal counts of 0 and 2) keeps the original value, matching the
    original `med == 1` branch.
    """
    n = len(s)
    if n == 0:
        return s
    end_dist = round(r / 10)

    s_arr = np.asarray(s, dtype=np.int32)

    # Prefix sum for O(1) window-sum queries: cum[i+1] - cum[i] = s_arr[i]
    cum = np.empty(n + 1, dtype=np.int64)
    cum[0] = 0
    np.cumsum(s_arr, out=cum[1:])

    idx = np.arange(n, dtype=np.int64)

    # Compute begin/end for every index in parallel, honouring the same
    # if-elif priority as the original.
    begins = np.where(idx < end_dist,       0,
             np.where(idx < r,              0,
             np.where(idx + end_dist >= n,  np.maximum(0, n - end_dist),
             np.where(idx + r >= n,         idx - (n - idx) + 1,
                                            idx - r))))
    ends   = np.where(idx < end_dist,       end_dist + 1,
             np.where(idx < r,              2 * idx + 1,
             np.where(idx + end_dist >= n,  n,
             np.where(idx + r >= n,         n,
                                            idx + r + 1))))

    begins = np.clip(begins, 0, n).astype(np.int64)
    ends   = np.clip(ends,   0, n).astype(np.int64)

    win_sums  = cum[ends] - cum[begins]
    win_sizes = ends - begins

    # majority 2 → 2, majority 0 → 0, tie → keep original
    s_new = np.where(win_sums > win_sizes, 2,
            np.where(win_sums < win_sizes, 0, s_arr))

    return s_new


def process_diffs(x, y, y_new):
    count = int(np.sum(np.asarray(y) != np.asarray(y_new)))
    print(f"\tpotentially wrong types: {count}", flush=True)


def break_points(x, y):
    y_arr = np.asarray(y)
    x_arr = np.asarray(x)
    break_ls = x_arr[1:][y_arr[1:] != y_arr[:-1]].tolist()
    print(f"\tbreakpoints: {len(break_ls)}", flush=True)
    return break_ls


def assess_blocks(x, y, truth_x, truth_y):
    """
    For each truth segment (defined by consecutive breakpoints in truth_x/y),
    find the test variants that fall within that segment and record whether
    each agrees (0) or disagrees (1) with the truth label.

    Uses np.searchsorted for O(log n) segment lookup instead of a Polars
    filter per breakpoint.  x must be sorted (genomic positions always are).
    """
    x_arr = np.asarray(x, dtype=np.int64)
    y_arr = np.asarray(y)

    new_x, new_y = [], []
    last_type  = truth_y[0]
    start_pos  = truth_x[0]

    for idx in range(1, len(truth_x)):
        r_x, r_y = truth_x[idx], truth_y[idx]
        if r_y != last_type:
            end_pos = truth_x[idx - 1]
            lo = np.searchsorted(x_arr, start_pos, side='left')
            hi = np.searchsorted(x_arr, end_pos,   side='right')
            seg_x = x_arr[lo:hi]
            seg_y = y_arr[lo:hi]
            new_x.extend(seg_x.tolist())
            new_y.extend((seg_y != last_type).astype(int).tolist())
            start_pos = r_x
        last_type = r_y

    lo = np.searchsorted(x_arr, start_pos,    side='left')
    hi = np.searchsorted(x_arr, truth_x[-1],  side='right')
    new_x.extend(x_arr[lo:hi].tolist())
    new_y.extend((y_arr[lo:hi] != last_type).astype(int).tolist())

    return new_x, new_y
