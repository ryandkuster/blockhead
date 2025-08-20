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


def deep_dive(args, format_fields, df, named_f1_dt):
    """
    Create new MIER column per-sample that is 1 (pass), 2 (fail), or 3
    (unknown).
    """
    print("analyzing parental trios with homozygous founder calls (deep dive mode)")
    dd_ls, dd_header = prepare_deep_dive_dt(format_fields)

    with open(os.path.join(args.outdir,"MIER_deep_dive_homozygous_parents.tsv"), "w") as o:
        o.write(f"{dd_header}\n")
        for f1, (p1, p2) in named_f1_dt.items():
            print(f"{f1} {p1} {p2}")
            tmp_df = df.select([f1, p1, p2])
            tmp_df = expand_dd_df(format_fields, tmp_df, [f1, p1, p2])

            new_col = f"MIER_{f1}"


            df_hom = tmp_df.filter(
                ((pl.col(f"{p1}_GT") == "0") & (pl.col(f"{p2}_GT") == "2")) |
                ((pl.col(f"{p2}_GT") == "0") & (pl.col(f"{p1}_GT") == "2"))
            )

            df_hom = df_hom.with_columns(
                pl.when((df_hom[f"{f1}_GT"].str.contains(r"\.")) | (df_hom[f"{p1}_GT"].str.contains(r"\.")) | (df_hom[f"{p2}_GT"].str.contains(r"\."))).then(3)
                .when((df_hom[f"{f1}_GT"] == "0") & (df_hom[f"{p1}_GT"] != "2") & (df_hom[f"{p2}_GT"] != "2")).then(1)
                .when((df_hom[f"{f1}_GT"] == "1") & ~((df_hom[f"{p1}_GT"] == "0") & (df_hom[f"{p2}_GT"] == "0")) & ~((df_hom[f"{p1}_GT"] == "2") & (df_hom[f"{p2}_GT"] == "2"))).then(1)
                .when((df_hom[f"{f1}_GT"] == "2") & (df_hom[f"{p1}_GT"] != "0") & (df_hom[f"{p2}_GT"] != "0")).then(1)
                .otherwise(2)
                .alias(new_col)
            )

            correct = df_hom.filter(pl.col(new_col) == 1).height
            incorrect = df_hom.filter(pl.col(new_col) == 2).height
            unknown = df_hom.filter(pl.col(new_col) == 3).height
            correct_known = correct/(correct+incorrect)
            correct_total = correct/(correct+incorrect+unknown)

            dd_avg_ls = process_averages(df_hom, dd_ls, f1, new_col)

            result = [f1, correct, incorrect, unknown, correct_known, correct_total] + dd_avg_ls
            result = [str(i) for i in result]
            result = '\t'.join(result)
            o.write(f"{result}\n")


def process_averages(df_hom, dd_ls, f1, new_col):
    dd_avg_ls = []
    for i in dd_ls:
        if i.endswith("_incorrect"):
            i = i[:-len("_incorrect")]
            i = f"{f1}_{i}"
            avg = df_hom.filter(pl.col(new_col) == 2).select(pl.col(i).str.to_integer().mean()).item()
        if i.endswith("_correct"):
            i = i[:-len("_correct")]
            i = f"{f1}_{i}"
            avg = df_hom.filter(pl.col(new_col) == 1).select(pl.col(i).str.to_integer().mean()).item()
        dd_avg_ls.append(avg)
    return dd_avg_ls


def expand_dd_df(format_fields, tmp_df, trio_ls):
    """
    Given a list of format fields, split the tmp_df trio format fields
    and recode on the fly
    """
    result_df = None

    for sample in trio_ls:
        for idx, i in enumerate(format_fields.split(":")):
            if i not in ["GT", "PL", "DP", "SP", "AD"]:
                continue

            header = f"{sample}_{i}"
            sample_df = tmp_df.select([
                pl.col(sample).str.split(":").list.get(idx).alias(header)
            ])

            if i == "GT":
                sample_df = recode_vcf(sample_df, [header])
                sample_df = recode_missing([header], sample_df)
            elif i == "PL":
                sample_df = sample_df.select(
                    pl.col(header).str.split(",").list.get(0).alias(f"{header}_0"),
                    pl.col(header).str.split(",").list.get(1).alias(f"{header}_1"),
                    pl.col(header).str.split(",").list.get(2).alias(f"{header}_2")
                )
            elif i == "AD":
                sample_df = sample_df.select(
                    pl.col(header).str.split(",").list.get(0).alias(f"{header}_0"),
                    pl.col(header).str.split(",").list.get(1).alias(f"{header}_1"),
                )

            # Append the sample_df to the growing result_df
            if result_df is None:
                result_df = sample_df
            else:
                result_df = result_df.with_columns(sample_df)
    return result_df


def prepare_deep_dive_dt(format_fields):
    dd_ls = []
    dd_header = ["f1", "correct", "incorrect", "unknown", "correct_known", "correct_total"]

    for i in format_fields.split(":"):
        if i == "GT":
            continue
        elif i == "PL":
            for j in range(0, 3):
                dd_ls.append(f"{i}_{j}_correct")
                dd_ls.append(f"{i}_{j}_incorrect")
        elif i == "AD":
            for j in range(0, 2):
                dd_ls.append(f"{i}_{j}_correct")
                dd_ls.append(f"{i}_{j}_incorrect")
        else:
            dd_ls.append(f"{i}_correct")
            dd_ls.append(f"{i}_incorrect")
    for i in dd_ls:
        dd_header.append(i)
    dd_header = "\t".join(dd_header)
    return dd_ls, dd_header


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


def wrongo_bongo(args, df_coords, df, named_f1_dt):
    """
    Create new MIER column per-sample that is coded using decimel
    representations of binary encoding _ bool traits:
        - 1: Member of trio is missing (1 is yes)
        - 2: Correct (1 is yes)
        - 4: F1 call is homozygous REF (0 can be het or homozygous alt)
        - 8: F1 call is homozygous ALT (0 can be het or homozygous ref)
        - 16: F1 call is heterozygous

    Examples:
        - "10010" (18) is not missing, correct, heterozygous
        - "00110" (6) is not missing, correct, homozygous reference
        - "01000" (8) is not missing, incorrect, homozygous alternate
        - "10001" (17) is missing, incorrect, heterozygous
    """
    print("analyzing parental trios for wrong calls")
    df = pl.concat([df_coords, df], how="horizontal")
    keep_ls = ["#CHROM", "POS"]
    mier_ls = []

    for f1, (p1, p2) in named_f1_dt.items():
        print(f"{f1} {p1} {p2}")
        new_col = f"MIER_{f1}"
        mier_ls.append(new_col)

        # Flip first bit if missing
        df = df.with_columns(
            pl.when((df[f1].str.contains(r"\.")) | (df[p1].str.contains(r"\.")) | (df[p2].str.contains(r"\."))).then(0b00001)
            .otherwise(0b00000)
            .alias(new_col)
        )

        # Flip second bit if correct
        df = df.with_columns(
            pl.when((df[f1].str.contains(r"\.")) | (df[p1].str.contains(r"\.")) | (df[p2].str.contains(r"\."))).then(pl.col(new_col) | 0b00000)
            .when((df[f1] == "0") & (df[p1] != "2") & (df[p2] != "2")).then(pl.col(new_col) | 0b00010)
            .when((df[f1] == "1") & ~((df[p1] == "0") & (df[p2] == "0")) & ~((df[p1] == "2") & (df[p2] == "2"))).then(pl.col(new_col) | 0b00010)
            .when((df[f1] == "2") & (df[p1] != "0") & (df[p2] != "0")).then(pl.col(new_col) | 0b00010)
            .otherwise(pl.col(new_col) | 0b00000)
            .alias(new_col)
        )

        # Flip third bit if homozygous reference
        df = df.with_columns(
            pl.when((df[f1] == "0")).then(pl.col(new_col) | 0b00100)
            .otherwise(pl.col(new_col) | 0b00000)
            .alias(new_col)
        )

        # Flip fourth bit if homozygous alternate
        df = df.with_columns(
            pl.when((df[f1] == "2")).then(pl.col(new_col) | 0b01000)
            .otherwise(pl.col(new_col) | 0b00000)
            .alias(new_col)
        )

        # Flip fifth bit if heterozygous
        df = df.with_columns(
            pl.when((df[f1] == "1")).then(pl.col(new_col) | 0b10000)
            .otherwise(pl.col(new_col) | 0b00000)
            .alias(new_col)
        )

    # These bitwise masks can be used to row-sum features.
    # First value is the must be 1 mask, second is must be zero mask
    # For example, incorrect, non-missing hets would be [16, 15]
    bitwise_dt = {"total":         [0, 0],
                  "missing":       [0b00001, 0],
                  "known":         [0, 0b00001],
                  "correct":       [0b00010, 0],
                  "incorrect":     [0, 0b00011],
                  "hom_ref":       [0b00100, 0],
                  "hom_alt":       [0b01000, 0],
                  "het":           [0b10000, 0],
                  "hom_ref_wrong": [0b00100, 0b11011],
                  "hom_alt_wrong": [0b01000, 0b10111],
                  "het_wrong":     [0b10000, 0b01111]}

    combination_expressions = []
    for mier_type, mask in bitwise_dt.items():
        keep_ls.append(mier_type)
        combination_expressions.append(
            pl.sum_horizontal([
                ((pl.col(c) & mask[0]) == mask[0]) & ((pl.col(c) & mask[1]) == 0)
                for c in mier_ls
            ]).alias(f"{mier_type}"))

    # Add the new columns to the DataFrame
    df = df.with_columns(
        *combination_expressions
    )

    # with pl.Config(tbl_cols=-1, tbl_rows=-1):
    # with pl.Config(tbl_cols=20, tbl_rows=60):
        # print(df)
    # sys.exit()

    keep_ls += mier_ls
    df = df.select(pl.col(keep_ls))

    return df, mier_ls, bitwise_dt


def plot_wrong_calls(args, mier_ls, df, bit_type, suffix):
    """
    Iterate one scaffold at a time to see the calls that are incorrect.
    """

    #alpha_val = 0.5
    alpha_val = 0.02
    chrom_ls = sorted(df["#CHROM"].unique().to_list())
    colors_dt = haplotype_colors(args)

    for chrom in chrom_ls:

        # filter to a single chromosome
        print(f"subsetting {chrom}")
        chrom_df = df.filter(
            (df["#CHROM"] == chrom)
        )

        # prepare plot outfile name and basis
        img_path = os.path.join(args.outdir, f"{chrom}_{len(mier_ls)}_wrong_calls_{suffix}.png")
        fig, axes = plt.subplots(nrows=len(mier_ls), ncols=2, sharex=True, figsize=(32, len(mier_ls)/3), gridspec_kw={'hspace': 0.3, 'width_ratios': [90, 1], 'wspace': 0.03})

        for idx, f1 in enumerate(mier_ls):
            # Create a plot.
            try:
                f1_color = colors_dt[f1[5:]]
            except KeyError:
                f1_color = "green"

            x = chrom_df["POS"].to_list()
            y = chrom_df[f1].to_list()
            new_x = []
            new_y = []
            colors_ls = []

            for x_val, y_val in zip(x, y):
                # confirm call is not missing and incorrect
                if y_val & 0b00011 == 0:
                    new_x.append(x_val)
                    new_y.append(1)
                    if y_val in bit_type:
                        colors_ls.append("red")
                    else:
                        colors_ls.append("grey")

            # axes[idx, 0].scatter(new_x, new_y, s=500, marker="|", c="black", edgecolor="none", alpha=alpha_val)
            axes[idx, 0].scatter(new_x, new_y, s=500, marker="|", c=colors_ls, edgecolor="none", alpha=alpha_val)
            axes[idx, 0].set_yticks([])

            axes[idx, 0].set_ylabel(f"{f1}", labelpad=100, va="center", ha="left", rotation=0)
            axes[idx, 1].set_facecolor(f1_color)

        # remove ticks and labels for the f1/p3 labels 
        for idx, i in enumerate(mier_ls):
            axes[idx, 1].set_yticklabels([])
            axes[idx, 1].set_yticks([])
            axes[idx, 1].set_xticklabels([])
            axes[idx, 1].set_xticks([])

        # remove the whitespace on the x axis that matplotlib defaults to
        for ax in axes:
            ax[0].autoscale(enable=True, axis='x', tight=True)

        axes[0, 0].set_title(f"{chrom}", pad=20)
        axes[0, 1].set_title("F1", pad=20)

        print(f"saving image to {img_path}")
        plt.savefig(img_path)


def plot_wrong_blocks(args, mier_ls, df, bitwise_dt):
    df = df.select(pl.exclude(mier_ls))

    df = calc_perecentages(df)

    features_dt = {"pct_hom_ref_wrong": ["incorrect 0/0 / total",   "#F6D55C"],
                   "pct_het_wrong":     ["incorrect 0/1 / total",   "#3CAEA3"],
                   "pct_hom_alt_wrong": ["incorrect 1/1 / total",   "#20639B"],
                   "pct_incorrect":     ["incorrect / total", "#ED553B"],
                   "pct_correct":       ["correct / total", "#173F5F"],
                   "pct_missing":       ["missing / total", "black"]}

    chrom_ls = sorted(df["#CHROM"].unique().to_list())

    # Output final df of window mode
    new_df = pl.DataFrame()

    for chrom in chrom_ls:

        new_chrom_df = pl.DataFrame()

        # filter to a single chromosome
        print(f"subsetting {chrom}")
        chrom_df = df.filter(
            (df["#CHROM"] == chrom)
        )

        max_pos = chrom_df.select(pl.col("POS").max()).item()
        x_ticks = [i for i in range(0, max_pos, 5000000)]
        tick_labels = [f'{x/1e6:.0f}Mb' if x != 0 else '0' for x in x_ticks]

        # prepare plot outfile name and basis
        img_path = os.path.join(args.outdir, f"{chrom}_{len(mier_ls)}_{str(args.breaks)}_VISUAL_SUMMARY.png")
        # fig, axes = plt.subplots(nrows=len(features_dt.keys()), ncols=1, sharex=True, sharey=True, figsize=(30, len(features_dt.keys())*3), gridspec_kw={'hspace': 0.2, 'wspace': 0.03})
        fig, axes = plt.subplots(nrows=len(features_dt.keys()), ncols=1, sharex=True, sharey=False, figsize=(30, len(features_dt.keys())*3), gridspec_kw={'hspace': 0.2, 'wspace': 0.03})

        for idx, (feature, val) in enumerate(features_dt.items()):
            averages, window_starts = block_windows(args, chrom_df, feature)
            if idx == 0:
                window_ends = [i+(args.breaks-1) for i in window_starts]
                chrom_label = [chrom for i in window_starts]
                new_chrom_df = new_chrom_df.with_columns(pl.Series("chrom", chrom_label))
                new_chrom_df = new_chrom_df.with_columns(pl.Series("start", window_starts))
                new_chrom_df = new_chrom_df.with_columns(pl.Series("end", window_ends))
            new_chrom_df = new_chrom_df.with_columns(pl.Series(feature, averages))
            axes[idx].plot(window_starts, averages, c=val[1])
            axes[idx].fill_between(window_starts, averages, color=val[1], alpha=0.7)
            axes[idx].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.2f}'))
            axes[idx].set_ylabel(f"{val[0]}", labelpad=150, va="center", ha="left", rotation=0, fontsize=14)
            axes[idx].margins(x=0.01, y=0.08)
        axes[0].set_title(f"{chrom}", pad=20)
        axes[0].set_xticks(x_ticks)
        axes[0].set_xticklabels(tick_labels)
        # axes[0].set_ylim(top=1)

        new_df = pl.concat([new_df, new_chrom_df])

        plt.tight_layout()
        plt.savefig(img_path)
    
    new_df.write_csv(os.path.join(args.outdir, f"{str(args.breaks)}_windows_wrong_calls.tsv"), separator="\t", include_header=True)


def calc_perecentages(df):
    df = df.with_columns(
        pl.when(pl.col("total") == 0)
        .then(0)
        .otherwise(pl.col("missing") / pl.col("total"))
        .alias("pct_missing")
    )
    df = df.with_columns(
        pl.when(pl.col("total") == 0)
        .then(0)
        .otherwise(pl.col("incorrect") / pl.col("total"))
        .alias("pct_incorrect")
    )
    df = df.with_columns(
        pl.when(pl.col("total") == 0)
        .then(0)
        .otherwise(pl.col("correct") / pl.col("total"))
        .alias("pct_correct")
    )
    df = df.with_columns(
        pl.when(pl.col("known") == 0)
        .then(0)
        .otherwise(pl.col("het_wrong") / pl.col("known"))
        .alias("pct_het_wrong")
    )
    df = df.with_columns(
        pl.when(pl.col("known") == 0)
        .then(0)
        .otherwise(pl.col("hom_ref_wrong") / pl.col("known"))
        .alias("pct_hom_ref_wrong")
    )
    df = df.with_columns(
        pl.when(pl.col("known") == 0)
        .then(0)
        .otherwise(pl.col("hom_alt_wrong") / pl.col("known"))
        .alias("pct_hom_alt_wrong")
    )
    return df


def block_windows(args, chrom_df, feature):
    pos_markers = chrom_df['POS'].to_numpy()
    feature_np = chrom_df[feature].to_numpy()

    # Define your window parameters
    window_size = args.breaks
    start_pos = 0
    max_pos = chrom_df.select(pl.col("POS").max()).item()
    end_pos = int(np.ceil(max_pos / window_size) * window_size)

    window_edges = np.arange(start_pos, end_pos + 1, window_size)
    window_starts = window_edges[:-1]
    window_ends = window_edges[1:] - 1

    averages = []
    for i in range(len(window_starts)):
        mask = (pos_markers >= window_starts[i]) & (pos_markers <= window_ends[i])
        if np.any(mask):
            avg = np.mean(feature_np[mask])
        else:
            avg = 0  # or np.nan if you prefer
        averages.append(avg)

    return averages, window_starts