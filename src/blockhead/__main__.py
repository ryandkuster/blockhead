# src/blockhead/__main__.py
import gzip
import os
import polars as pl
import sys

import blockhead.df_manip as dm
import blockhead.parentage as pm
import blockhead.utils as um

def main():
    args = um.parse_user_input()
    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    og_df, df, df_coords = dm.read_vcf(args)
    parent_dt, cross_ls, named_f1_dt = pm.get_parentage(args)
    adv_ls = pm.get_advanced(named_f1_dt)
    sample_ls = pm.get_sample_ls(named_f1_dt)
    df = dm.recode_vcf(df, sample_ls)
    df = dm.recode_missing(sample_ls, df)

    df, mier_ls = dm.parental_trios(args, sample_ls, df, named_f1_dt)
    df = df.with_columns(
        (
            pl.sum_horizontal([pl.col(c) == 1 for c in mier_ls])
            / len(mier_ls)
        ).alias("percent_mier_correct")
    )
    df = dm.stitch_coords(df, df_coords)
    df = df.filter(df["percent_mier_correct"] >= args.threshold)

    if args.outvcf:
        df_coords = df.select(["#CHROM", "POS"])
        compressed = um.gzip_test(args.vcf)
        if compressed:
            with gzip.open(args.vcf, "rt", encoding='utf-8') as f:
                dm.write_outfile(args, og_df, df_coords, f)
        else:
            with open(args.vcf) as f:
                dm.write_outfile(args, og_df, df_coords, f)

    if args.blockmode and len(adv_ls) > 0:
        print("determining parentage haplotypes for advance hybrids")
        dm.haplotype_per_cross(args, adv_ls, named_f1_dt, df)


if __name__ == "__main__":
    main()
