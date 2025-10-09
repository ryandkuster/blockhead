import os

import matplotlib.pyplot as plt
import numpy as np
import polars as pl


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

        # Flip second bit if correct, if missing, don't modify this bit.
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
    df = df.with_columns(
        pl.col("#CHROM").str.split("#").list.get(-1).alias("#CHROM")
    )
    df.write_csv(os.path.join(args.outdir, f"MIER_per_snp_summary.tsv"), separator="\t", include_header=True)

    return df, mier_ls, bitwise_dt


def plot_wrong_blocks(args, mier_ls, df, bitwise_dt):
    """
    Create blocks (windows)
    """
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

        # Filter to a single chromosome
        print(f"subsetting {chrom}")
        chrom_df = df.filter(
            (df["#CHROM"] == chrom)
        )

        max_pos = chrom_df.select(pl.col("POS").max()).item()
        x_ticks = [i for i in range(0, max_pos, 5000000)]
        tick_labels = [f'{x/1e6:.0f}Mb' if x != 0 else '0' for x in x_ticks]

        # Prepare plot outfile name and basis
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
    """
    For a given feature, calculate averages within the defined block.
    """
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
            avg = 0
        averages.append(avg)

    return averages, window_starts