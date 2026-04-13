import os

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from blockhead.df_manip import HOM_REF, HET, UNKNOWN, HOM_ALT


def compute_wrong_calls_chrom(chrom, pos_arr, gts_matrix, sample_idx_in_matrix,
                               named_f1_dt):
    """
    Compute per-SNP bitwise MIER summary for one chromosome.
    Returns a Polars DataFrame with #CHROM, POS, summary columns, and MIER_* columns.

    Bit encoding:
        bit 1 (0b00001): missing
        bit 2 (0b00010): correct (and not missing)
        bit 3 (0b00100): F1 is hom_ref
        bit 4 (0b01000): F1 is hom_alt
        bit 5 (0b10000): F1 is het
    """
    n = len(pos_arr)
    mier_cols = {}

    for f1, (p1, p2) in named_f1_dt.items():
        f1_col = gts_matrix[:, sample_idx_in_matrix[f1]]
        p1_col = gts_matrix[:, sample_idx_in_matrix[p1]]
        p2_col = gts_matrix[:, sample_idx_in_matrix[p2]]

        missing = ((f1_col == UNKNOWN) | (p1_col == UNKNOWN) | (p2_col == UNKNOWN))

        correct = (
            ((f1_col == HOM_REF) & (p1_col != HOM_ALT) & (p2_col != HOM_ALT)) |
            ((f1_col == HET) &
             ~((p1_col == HOM_REF) & (p2_col == HOM_REF)) &
             ~((p1_col == HOM_ALT) & (p2_col == HOM_ALT))) |
            ((f1_col == HOM_ALT) & (p1_col != HOM_REF) & (p2_col != HOM_REF))
        ) & ~missing

        bitwise = (
            missing.astype(np.uint8)             * 0b00001 +
            correct.astype(np.uint8)             * 0b00010 +
            (f1_col == HOM_REF).astype(np.uint8) * 0b00100 +
            (f1_col == HOM_ALT).astype(np.uint8) * 0b01000 +
            (f1_col == HET).astype(np.uint8)     * 0b10000
        )
        mier_cols[f"MIER_{f1}"] = bitwise.tolist()

    mier_ls = list(mier_cols.keys())

    df_dict = {"#CHROM": [chrom] * n, "POS": pos_arr.tolist()}
    df_dict.update(mier_cols)
    df = pl.DataFrame(df_dict)

    bitwise_dt = {
        "total":         [0,        0],
        "missing":       [0b00001,  0],
        "known":         [0,        0b00001],
        "correct":       [0b00010,  0],
        "incorrect":     [0,        0b00011],
        "hom_ref":       [0b00100,  0],
        "hom_alt":       [0b01000,  0],
        "het":           [0b10000,  0],
        "hom_ref_wrong": [0b00100,  0b11011],
        "hom_alt_wrong": [0b01000,  0b10111],
        "het_wrong":     [0b10000,  0b01111],
    }

    exprs = []
    for mier_type, (must_one, must_zero) in bitwise_dt.items():
        exprs.append(
            pl.sum_horizontal([
                ((pl.col(c) & must_one) == must_one) & ((pl.col(c) & must_zero) == 0)
                for c in mier_ls
            ]).alias(mier_type)
        )
    df = df.with_columns(*exprs)

    keep = (["#CHROM", "POS"] + list(bitwise_dt.keys()) + mier_ls)
    return df.select(keep)


def plot_wrong_blocks(args, mier_ls, df):
    df = df.select(pl.exclude(mier_ls))
    df = calc_perecentages(df)

    features_dt = {
        "pct_hom_ref_wrong": ["incorrect 0/0 / total",  "#F6D55C"],
        "pct_het_wrong":     ["incorrect 0/1 / total",  "#3CAEA3"],
        "pct_hom_alt_wrong": ["incorrect 1/1 / total",  "#20639B"],
        "pct_incorrect":     ["incorrect / total",      "#ED553B"],
        "pct_correct":       ["correct / total",        "#173F5F"],
        "pct_missing":       ["missing / total",        "black"],
    }

    chrom_ls = sorted(df["#CHROM"].unique().to_list())
    new_df = pl.DataFrame()

    for chrom in chrom_ls:
        print(f"subsetting {chrom}", flush=True)
        chrom_df = df.filter(pl.col("#CHROM") == chrom)
        max_pos = chrom_df.select(pl.col("POS").max()).item()
        x_ticks = list(range(0, max_pos, 5_000_000))
        tick_labels = [f'{x/1e6:.0f}Mb' if x != 0 else '0' for x in x_ticks]

        img_path = os.path.join(
            args.outdir,
            f"{chrom}_{len(mier_ls)}_{str(args.breaks)}_wrong_calls.png"
        )
        fig, axes = plt.subplots(
            nrows=len(features_dt) + 1, ncols=1, sharex=True, sharey=False,
            figsize=(30, (len(features_dt) + 1) * 3),
            gridspec_kw={'hspace': 0.2, 'wspace': 0.03}
        )

        new_chrom_df = pl.DataFrame()
        for idx, (feature, (label, color)) in enumerate(features_dt.items()):
            averages, counts, window_starts = block_windows(args, chrom_df, feature)
            if idx == 0:
                window_ends = [i + (args.breaks - 1) for i in window_starts]
                new_chrom_df = new_chrom_df.with_columns(
                    pl.Series("chrom", [chrom] * len(window_starts)),
                    pl.Series("start", window_starts),
                    pl.Series("end", window_ends),
                )
            new_chrom_df = new_chrom_df.with_columns(pl.Series(feature, averages))
            axes[idx].plot(window_starts, averages, c=color)
            axes[idx].fill_between(window_starts, averages, color=color, alpha=0.7)
            axes[idx].yaxis.set_major_formatter(
                plt.FuncFormatter(lambda x, p: f'{x:.2f}'))
            axes[idx].set_ylabel(label, labelpad=150, va="center",
                                 ha="left", rotation=0, fontsize=14)
            axes[idx].margins(x=0.01, y=0.08)

        axes[idx + 1].plot(window_starts, counts, c="grey")
        axes[idx + 1].fill_between(window_starts, counts, color="grey", alpha=0.7)
        axes[idx + 1].yaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, p: f'{x:.2f}'))
        axes[idx + 1].set_ylabel("SNP count", labelpad=150, va="center",
                                  ha="left", rotation=0, fontsize=14)
        axes[idx + 1].margins(x=0.01, y=0.08)

        axes[0].set_title(f"{chrom}", pad=20)
        axes[0].set_xticks(x_ticks)
        axes[0].set_xticklabels(tick_labels)

        new_df = pl.concat([new_df, new_chrom_df])
        plt.tight_layout()
        plt.savefig(img_path)
        plt.close(fig)

    new_df.write_csv(
        os.path.join(args.outdir, f"{str(args.breaks)}_windows_wrong_calls.tsv"),
        separator="\t", include_header=True
    )


def calc_perecentages(df):
    df = df.with_columns(
        pl.when(pl.col("total") == 0).then(0)
          .otherwise(pl.col("missing") / pl.col("total"))
          .alias("pct_missing"),
        pl.when(pl.col("total") == 0).then(0)
          .otherwise(pl.col("incorrect") / pl.col("total"))
          .alias("pct_incorrect"),
        pl.when(pl.col("total") == 0).then(0)
          .otherwise(pl.col("correct") / pl.col("total"))
          .alias("pct_correct"),
        pl.when(pl.col("known") == 0).then(0)
          .otherwise(pl.col("het_wrong") / pl.col("total"))
          .alias("pct_het_wrong"),
        pl.when(pl.col("known") == 0).then(0)
          .otherwise(pl.col("hom_ref_wrong") / pl.col("total"))
          .alias("pct_hom_ref_wrong"),
        pl.when(pl.col("known") == 0).then(0)
          .otherwise(pl.col("hom_alt_wrong") / pl.col("total"))
          .alias("pct_hom_alt_wrong"),
    )
    return df


def block_windows(args, chrom_df, feature):
    pos_markers = chrom_df['POS'].to_numpy()
    feature_np = chrom_df[feature].to_numpy()

    window_size = args.breaks
    max_pos = chrom_df.select(pl.col("POS").max()).item()
    end_pos = int(np.ceil(max_pos / window_size) * window_size)

    window_edges = np.arange(0, end_pos + 1, window_size)
    window_starts = window_edges[:-1]
    window_ends = window_edges[1:] - 1

    averages, counts = [], []
    for ws, we in zip(window_starts, window_ends):
        mask = (pos_markers >= ws) & (pos_markers <= we)
        if np.any(mask):
            averages.append(float(np.mean(feature_np[mask])))
            counts.append(int(np.sum(mask)))
        else:
            averages.append(0.0)
            counts.append(0)

    return averages, counts, window_starts.tolist()
