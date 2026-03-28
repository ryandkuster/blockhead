#!/usr/bin/env python3

"""
The inputs are:
    sys.argv[1] is linear summary of windows (wrong calls) from blockhead
    sys.argv[2] is graph summary of windows (wrong calls) from blockhead
    sys.argv[3] is the output directory
"""
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
import polars as pl
import sys


def main():
    features_dt = {"pct_hom_ref_wrong": [f"Incorrect\n(0/0)", "#F6D55C"],
                   "pct_het_wrong":     [f"Incorrect\n(0/1)", "#3CAEA3"],
                   "pct_hom_alt_wrong": [f"Incorrect\n(1/1)", "#20639B"],
                   "pct_incorrect":     [f"Incorrect\n(all)", "#ED553B"],
                   "pct_correct":       [f"Correct", "#173F5F"],
                   "pct_missing":       [f"Missing", "black"]}

    df1 = pl.read_csv(sys.argv[1], separator="\t")
    df2 = pl.read_csv(sys.argv[2], separator="\t")
    print(df1)
    print(df2)

    chrom_ls1 = sorted(df1["chrom"].unique().to_list())
    chrom_ls2 = sorted(df2["chrom"].unique().to_list())

    for chrom1, chrom2 in zip(chrom_ls1, chrom_ls2):
        # filter to a single chromosome
        print(f"subsetting {chrom1} and {chrom2}")
        chrom_df1 = df1.filter(
            (df1["chrom"] == chrom1)
        )
        chrom_df2 = df2.filter(
            (df2["chrom"] == chrom2)
        )

        img_path = os.path.join(sys.argv[3], f"compare_{chrom1}.png")
        fig, axes = plt.subplots(nrows=len(features_dt.keys()), ncols=1, sharex=True, sharey=False, figsize=(30, len(features_dt.keys())*3), gridspec_kw={'hspace': 0.2, 'wspace': 0.03})

        for idx, (feature, val) in enumerate(features_dt.items()):
            axes[idx].plot(chrom_df1["start"], chrom_df1[feature] * 100, c=val[1], linewidth=3)
            axes[idx].plot(chrom_df2["start"], chrom_df2[feature] * 100, c=val[1], linewidth=3, linestyle='dotted')
            axes[idx].fill_between(chrom_df1["start"], chrom_df1[feature] * 100, chrom_df2[feature] * 100, alpha=0.1, color="gray")
            axes[idx].set_ylabel(f"{val[0]}", va="center", ha="center", rotation=0, fontsize=22)
            axes[idx].yaxis.set_label_coords(-0.08, 0.5)
            axes[idx].margins(x=0.01, y=0.08)
            axes[idx].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.0f}%'))
            axes[idx].tick_params(axis='y', labelsize=14)
            axes[idx].tick_params(axis='x', length=5, width=2)

        max_pos = chrom_df1.select(pl.col("start").max()).item()
        x_ticks = [i for i in range(0, max_pos, 5000000)]
        tick_labels = [f'{x/1e6:.0f}Mb' if x != 0 else '0' for x in x_ticks]
        axes[0].set_title(f"{chrom1}", pad=20, fontsize=30)
        axes[-1].set_xticks(x_ticks)
        axes[-1].set_xticklabels(tick_labels)
        axes[-1].tick_params(axis='x', labelsize=18)

        graph = mlines.Line2D([], [], color='black', linestyle='dotted', label='vg-surject BCFtools', linewidth=3)
        linear = mlines.Line2D([], [], color='black', linestyle='solid', label='linear BCFtools', linewidth=3)

        plt.legend(
            bbox_to_anchor=(1, len(features_dt.keys())+1),
            loc='upper right',
            frameon=True,       # show frame
            shadow=True,        # add shadow
            fancybox=True,      # rounded corners
            ncol=2,             # number of columns
            handles=[graph, linear],
            fontsize=22
        )
        plt.tight_layout()
        plt.savefig(img_path, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
