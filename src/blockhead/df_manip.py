import gzip
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns

from blockhead.parentage import get_advanced_lineage

# pl.Config.set_tbl_cols(-1)

def read_vcf(args) -> pl.DataFrame:
    """
    Open the vcf, skipping all header info except relevant fields.
    Keep only bi-allelic SNPs.
    Save the df_coords for writing final output file.
    """
    print("opening vcf")
    og_df = pl.read_csv(args.vcf,
                     separator="\t",
                     comment_prefix="##")
    print(f"{og_df.shape[0]} variants found")
    df = og_df.filter((pl.col("REF").str.len_chars() == 1) & (pl.col("ALT").str.len_chars() == 1))
    print(f"{df.shape[0]} biallelic SNPs found")
    df_coords = df.select(["#CHROM", "POS"])
    format_fields = df.select(pl.col(df.columns[8]).first()).item(0, 0)
    df = df.select(df.columns[9:])
    return og_df, df, df_coords, format_fields


def recode_vcf(df: pl.DataFrame, sample_ls: dict) -> pl.DataFrame:
    """
    Recode three monoallelic genotypes as 0, 1, or 2.
    """
    df = df.with_columns([
        pl.col(sample).str.split(":").list.first()
        .str.replace(r"[\/|]", "", literal=False)
        .str.replace(r"00", "0")
        .str.replace(r"11", "2")
        .str.replace(r"01", "1")
        .str.replace(r"10", "1").alias(sample).cast(pl.Int64, strict=False)
    for sample in sample_ls
    ])
    return df


def recode_missing(sample_ls, df):
    """
    If multiallelic sites present, recode as ".." unknown
    """
    df = df.with_columns(
        pl.col(sample_ls).fill_null("..")
    )
    return df


def parental_trios(args, sample_ls, df, named_f1_dt):
    """
    Create new MIER column per-sample that is 1 (pass), 2 (fail), or 3
    (unknown).
    """
    print("analyzing parental trios")
    mier_ls = []
    
    with open(os.path.join(args.outdir,"MIER_summary.tsv"), "w") as o:
        o.write(f"f1\tcorrect\tincorrect\tunknown\tcorrect_known\tcorrect_total\n")
        for f1, (p1, p2) in named_f1_dt.items():
            print(f"{f1} {p1} {p2}")
            new_col = f"MIER_{f1}"
            mier_ls.append(new_col)
            sample_ls.append(new_col)

            df = df.with_columns(
                pl.when((df[f1].str.contains(r"\.")) | (df[p1].str.contains(r"\.")) | (df[p2].str.contains(r"\."))).then(3)
                .when((df[f1] == "0") & (df[p1] != "2") & (df[p2] != "2")).then(1)
                .when((df[f1] == "1") & ~((df[p1] == "0") & (df[p2] == "0")) & ~((df[p1] == "2") & (df[p2] == "2"))).then(1)
                .when((df[f1] == "2") & (df[p1] != "0") & (df[p2] != "0")).then(1)
                .otherwise(2)
                .alias(new_col)
            )
            correct = df.filter(pl.col(new_col) == 1).height
            incorrect = df.filter(pl.col(new_col) == 2).height
            unknown = df.filter(pl.col(new_col) == 3).height

            o.write(f"{f1}\t{correct}\t{incorrect}\t{unknown}\t{(correct/(correct+incorrect)):.3f}\t{(correct/(correct+incorrect+unknown)):.3f}\n")

    return df, mier_ls


def homozygous_parents(args, df, named_f1_dt):
    """
    Create new MIER column per-sample that is 1 (pass), 2 (fail), or 3
    (unknown).
    """
    print("analyzing parental trios with homozygous founder calls")

    with open(os.path.join(args.outdir,"MIER_summary_homozygous_parents.tsv"), "w") as o:
        o.write(f"f1\tcorrect\tincorrect\tunknown\tcorrect_known\tcorrect_total\n")
        for f1, (p1, p2) in named_f1_dt.items():
            print(f"{f1} {p1} {p2}")
            new_col = f"MIER_{f1}"

            df_hom = df.filter(
                ((pl.col(p1) == "0") & (pl.col(p2) == "2")) |
                ((pl.col(p2) == "0") & (pl.col(p1) == "2"))
            )

            df_hom = df_hom.with_columns(
                pl.when((df_hom[f1].str.contains(r"\.")) | (df_hom[p1].str.contains(r"\.")) | (df_hom[p2].str.contains(r"\."))).then(3)
                .when((df_hom[f1] == "0") & (df_hom[p1] != "2") & (df_hom[p2] != "2")).then(1)
                .when((df_hom[f1] == "1") & ~((df_hom[p1] == "0") & (df_hom[p2] == "0")) & ~((df_hom[p1] == "2") & (df_hom[p2] == "2"))).then(1)
                .when((df_hom[f1] == "2") & (df_hom[p1] != "0") & (df_hom[p2] != "0")).then(1)
                .otherwise(2)
                .alias(new_col)
            )
            correct = df_hom.filter(pl.col(new_col) == 1).height
            incorrect = df_hom.filter(pl.col(new_col) == 2).height
            unknown = df_hom.filter(pl.col(new_col) == 3).height

            o.write(f"{f1}\t{correct}\t{incorrect}\t{unknown}\t{(correct/(correct+incorrect)):.3f}\t{(correct/(correct+incorrect+unknown)):.3f}\n")


def stitch_coords(df, df_coords):
    df = pl.concat([df_coords, df], how="horizontal")
    return df


def write_outfile(args, og_df, df_coords, f):
    """
    Write a vcf of the input vcf file coords with only the relevant
    MIER correct variants present based on filtering threshold.
    """

    original_count = og_df.shape[0]
    og_df = og_df.join(df_coords, on=["#CHROM", "POS"], how="semi")
    percent_count = og_df.shape[0]/original_count
    args.outfile = os.path.join(args.outdir, "MIER_filtered.vcf")

    print("writing header")
    with open(args.outfile, "w") as o:
        for line in f:
            if line.startswith("##"):
                o.write(line)
            else:
                break

    print(f"writing {og_df.shape[0]} variants ({percent_count:.2%} original vcf entries retained)")
    with open(args.outfile, "a") as o:
        og_df.write_csv(o, separator='\t')


def haplotype_per_cross(args, adv_ls, named_f1_dt, df):
    """
    Iterate one scaffold at a time to see the calls that likely were
    contributed by each original parental cross (p1/p2).

    - [x] filter to p1/p2 == 0/2 or 2/0
    - [x] filter to f1 mier correct
    - [x] filter to p3 homozygous
    """

    # alpha_val = 0.01
    alpha_val = 0.5
    chrom_ls = sorted(df["#CHROM"].unique().to_list())
    colors_dt = haplotype_colors(args)

    for chrom in chrom_ls:

        # filter to a single chromosome
        print(f"subsetting {chrom}")
        chrom_df = df.filter(
            (df["#CHROM"] == chrom)
        )
        img_path = os.path.join(args.outdir, f"{chrom}_{len(adv_ls)}_haplotypes.png")
        fig, axes = plt.subplots(nrows=len(adv_ls), ncols=3, sharex=True, figsize=(32, len(adv_ls)/3), gridspec_kw={'hspace': 0.3, 'width_ratios': [90, 1, 1], 'wspace': 0.03})


        for idx, adv in enumerate(adv_ls):
            print(f"processing advanced hybrid : {adv}")
            lin_dt = get_advanced_lineage(adv, named_f1_dt)
            print(lin_dt)

            # filter to p1/p2 0/2 or 2/0
            tmp_df = chrom_df.filter(
                ((chrom_df[lin_dt["p1"]] == "0") & (chrom_df[lin_dt["p2"]] == "2")) | \
                ((chrom_df[lin_dt["p1"]] == "2") & (chrom_df[lin_dt["p2"]] == "0"))
            )

            # filter to f1 MIER correct
            tmp_df = tmp_df.filter(
                (tmp_df[f"MIER_{lin_dt["f1"]}"] == "1")
            )

            # filter to f1 MIER correct
            tmp_df = tmp_df.filter(
                (tmp_df[f"MIER_{lin_dt["adv"]}"] == "1")
            )

            # filter to p3 homozygous
            tmp_df = tmp_df.filter(
                (tmp_df[lin_dt["p3"]] == "0") | \
                (tmp_df[lin_dt["p3"]] == "2")
            )

            """
            Create new binary column called parent_type.
            adv p3  parent_type
            --- --  ------
            0   0   0
            1   0   2
            1   2   0
            2   2   2
            """
            tmp_df = tmp_df.with_columns(
                pl.when((tmp_df[lin_dt["adv"]] == "0") & (tmp_df[lin_dt["p3"]] == "0")).then(pl.lit("0"))
                  .when((tmp_df[lin_dt["adv"]] == "1") & (tmp_df[lin_dt["p3"]] == "0")).then(pl.lit("2"))
                  .when((tmp_df[lin_dt["adv"]] == "1") & (tmp_df[lin_dt["p3"]] == "2")).then(pl.lit("0"))
                  .when((tmp_df[lin_dt["adv"]] == "2") & (tmp_df[lin_dt["p3"]] == "2")).then(pl.lit("2"))
                  .otherwise(pl.lit("99"))
                  .alias("parent_type")
            )

            """
            Create new column called parent_origin. This connects the
            parent_type with the origin parent.
            """
            tmp_df = tmp_df.with_columns(
                pl.when((tmp_df["parent_type"] == tmp_df[lin_dt["p1"]])).then(pl.lit(lin_dt["p1"]))
                  .when((tmp_df["parent_type"] == tmp_df[lin_dt["p2"]])).then(pl.lit(lin_dt["p2"]))
                  .otherwise(99)
                  .alias("parent_origin")
            )

            keep_cols = ["POS", "parent_origin"] + [v for k,v in lin_dt.items()]
            tmp_df = tmp_df.select(keep_cols)

            tmp_df = tmp_df.with_columns(
                pl.col("parent_origin").cast(pl.Categorical).to_physical().alias("parent_cats")
            )

            # Create a plot.
            x = tmp_df["POS"].to_list()
            y = tmp_df["parent_origin"].to_list()
            colors = [colors_dt[id] for id in y]
            y = [0 for i in y]
            # axes[idx].set_ylabel(f"{adv}", labelpad=40, loc="center", rotation=0)
            axes[idx, 0].scatter(x, y, s=500, marker="|", c=colors, edgecolor="none", alpha=alpha_val)
            axes[idx, 0].set_yticks([])
            # axes[idx, 0].set_ylabel(f"{adv}", labelpad=50, loc="center", rotation=0)
            axes[idx, 0].set_ylabel(f"{adv}", labelpad=100, va="center", ha="left", rotation=0)

            axes[idx, 1].set_facecolor(colors_dt[lin_dt["f1"]])
            axes[idx, 2].set_facecolor(colors_dt[lin_dt["p3"]])

        # remove ticks and labels for the f1/p3 labels 
        for idx, i in enumerate(adv_ls):
            axes[idx, 1].set_yticklabels([])
            axes[idx, 1].set_yticks([])
            axes[idx, 1].set_xticklabels([])
            axes[idx, 1].set_xticks([])
            axes[idx, 2].set_yticklabels([])
            axes[idx, 2].set_yticks([])
            axes[idx, 2].set_xticklabels([])
            axes[idx, 2].set_xticks([])

        # remove the whitespace on the x axis that matplotlib defaults to
        for ax in axes:
            ax[0].autoscale(enable=True, axis='x', tight=True)

        axes[0, 0].set_title(f"{chrom}", pad=20)
        axes[0, 1].set_title("F1", pad=20)
        axes[0, 2].set_title("P3", pad=20)

        print(f"saving image to {img_path}")
        plt.savefig(img_path)

    return tmp_df


def haplotype_colors(args):
    if args.colors:
        colors_df = pl.read_csv(args.colors, separator="\t", has_header=False, new_columns=["id", "color"])
        colors_dt = dict(zip(colors_df["id"], colors_df["color"]))
    else:
        colors_dt = {}

    return colors_dt
