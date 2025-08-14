import gzip
import os
import sys

import matplotlib.pyplot as plt
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
        #fig, axes = plt.subplots(nrows=len(adv_ls), sharex=True, figsize=(30, len(adv_ls)/3), gridspec_kw={'hspace': 0.3})
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
    representations of binary encoding 4 bool traits:
        - Both alleles wrong relative to parents (1 means both, 0 is 1)
        - F1 call is homozygous REF (0 can be het or homozygous alt)
        - F1 call is homozygous ALT (0 can be het or homozygous ref)
        - Call is incorrect in the direction of REF
    
    Example "1010" (represented in df as 10):
        - Wrong relative to parents at both alleles
        - Homozygous REF is False
        - Homozygous ALT is True
        - Not biased in direction of REF (F1 is 1/1; parents are both 0/0)

    For all calls that are MIER correct or have a missing call for a
    trio member, codea as -1
    """
    print("analyzing parental trios for wrong calls")
    df = pl.concat([df_coords, df], how="horizontal")
    keep_ls = ["#CHROM", "POS"]
    mier_ls = []

    for f1, (p1, p2) in named_f1_dt.items():
        print(f"{f1} {p1} {p2}")
        new_col = f"MIER_{f1}"
        # new_col = f1
        keep_ls.append(new_col)
        mier_ls.append(new_col)

        # using (number_wrong)(hom_ref)(hom_alt)(ref_biased)
        df = df.with_columns(
            pl.when((df[f1].str.contains(r"\.")) | (df[p1].str.contains(r"\.")) | (df[p2].str.contains(r"\."))).then(-1)
            .when((df[f1] == "0") & (df[p1] != "2") & (df[p2] != "2")).then(-1)
            .when((df[f1] == "1") & ~((df[p1] == "0") & (df[p2] == "0")) & ~((df[p1] == "2") & (df[p2] == "2"))).then(-1)
            .when((df[f1] == "2") & (df[p1] != "0") & (df[p2] != "0")).then(-1)
            .when((df[p1] == "0") & (df[p2] == "0") & (df[f1] == "1")).then(0)  # 0000 0 
            .when((df[p1] == "0") & (df[p2] == "0") & (df[f1] == "2")).then(10) # 1010 10
            .when((df[p1] == "0") & (df[p2] == "1") & (df[f1] == "2")).then(2)  # 0010 2
            .when((df[p1] == "0") & (df[p2] == "2") & (df[f1] == "0")).then(5)  # 0101 5
            .when((df[p1] == "0") & (df[p2] == "2") & (df[f1] == "2")).then(2)  # 0010 2
            .when((df[p1] == "1") & (df[p2] == "0") & (df[f1] == "2")).then(2)  # 0010 2
            .when((df[p1] == "1") & (df[p2] == "2") & (df[f1] == "0")).then(5)  # 0101 5 
            .when((df[p1] == "2") & (df[p2] == "0") & (df[f1] == "0")).then(5)  # 0101 5 
            .when((df[p1] == "2") & (df[p2] == "0") & (df[f1] == "2")).then(2)  # 0010 2
            .when((df[p1] == "2") & (df[p2] == "1") & (df[f1] == "0")).then(5)  # 0101 5
            .when((df[p1] == "2") & (df[p2] == "2") & (df[f1] == "0")).then(13) # 1101 13
            .when((df[p1] == "2") & (df[p2] == "2") & (df[f1] == "1")).then(1)  # 0001 1
            .otherwise(999)
            .alias(new_col)
        )

    df = df.select(pl.col(keep_ls))

    # remove calls where all samples are correct or missing
    df = df.filter(
        ~pl.all_horizontal([pl.col(c) < 0 for c in mier_ls])
    )

    df = df.with_columns(
        pl.sum_horizontal([pl.col(c) >= 0 for c in mier_ls])
        .alias("sum_incorrect")
    )

    return df, mier_ls


def plot_wrong_calls(args, mier_ls, df, type_ls, suffix):
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
        img_path = os.path.join(args.outdir, f"{chrom}_{len(mier_ls)}_wrong_calls_{suffix}.png")
        fig, axes = plt.subplots(nrows=len(mier_ls), ncols=2, sharex=True, figsize=(32, len(mier_ls)/3), gridspec_kw={'hspace': 0.3, 'width_ratios': [90, 1], 'wspace': 0.03})

        for idx, f1 in enumerate(mier_ls):
            print(f"processing : {f1}")

            # Create a plot.
            f1_color = colors_dt[f1[5:]]
            x = chrom_df["POS"].to_list()
            y = chrom_df[f1].to_list()
            new_x = []
            new_y = []
            colors_ls = []

            for x_val, y_val in zip(x, y):
                if y_val >= 0:
                    new_x.append(x_val)
                    new_y.append(1)
                    # if y_val in [10, 13]: # wrong both alleles
                    # if y_val in [5, 13]: # homozygous REF
                    # if y_val in [2, 10]: # homozygous ALT
                    if y_val in type_ls:
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
