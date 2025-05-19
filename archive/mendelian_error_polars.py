#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "polars",
# ]
# ///

"""
Notes:
    - If a parent is het, calculate inheritance rate in F1s (threshold for
      removal if not at all near .5).
    - The order of 0/1 is arbitrary for non-phased but maybe write it as though
      phasing could come later. By using a window approach we can determine
      inheritance to be off if nearby variants are assorting at a different
      rate into F1.

Inputs:
    - vcf file with parents and F1s
    - file showing the expected combos of parents to produce F1s (3 columns)

Confirm f1 genotype is possible given parental genotypes:
    - 0/0 cannot originate from p1 = 1/1 OR p2 = 1/1
    - 0/1 cannot originate from (p1 = 1/1 and p2 = 1/1) OR (p1 = 0/0 and p2 = 0/0)
    - 1/1 cannot originate from p1 = 0/0 OR p2 = 0/0
    - order of 0/1 vs. 1/0 is not treated as phased
    - multiallelic genotypes (alleles = 2+) or missing (./1 etc.) in f1, p1, or p2 counted as unknown
"""

import gzip
import polars as pl
import re
import sys

def main():
    parentage = sys.argv[1] # parentage file
    vcf = sys.argv[2]       # gzipped vcf file
    results = sys.argv[3]   # tsv summary file

    if len(sys.argv) > 4:
        outfile = sys.argv[4]
    else:
        outfile = None

    named_f1_dt = get_parentage(parentage)

    sample_ls = []
    for k, v in named_f1_dt.items():
        sample_ls.append(k)
        for i in v:
            sample_ls.append(i)
    sample_ls = list(set(sample_ls))

    print("opening vcf")
    df = pl.scan_csv(vcf,
                     separator="\t",
                     comment_prefix="##")
    df = df.collect()
    print(df)
    print("removing multi-allelic sites, if present")
    df = df.filter(~pl.col("ALT").str.contains(","))
    print(df)
    df_coords = df.select(["#CHROM", "POS"])
    df = df.select(df.columns[9:])

    for sample in sample_ls:

        # recode three monoallelic genotypes as 0, 1, or 2
        print("recoding variants")
        df = df.with_columns(
            pl.col(sample).str.split(":").list.first()
            .str.replace(r"[\/|]", "", literal=False)
            .str.replace(r"00", "0")
            .str.replace(r"11", "2")
            .str.replace(r"01", "1")
            .str.replace(r"10", "1").alias(sample)
        )

    print("analyzing parental trios")
    with open(results, "w") as o:
        o.write(f"f1\tcorrect\tincorrect\tunknown\tcorrect_known\tcorrect_total\n")
        for f1, (p1, p2) in named_f1_dt.items():
            new_col = f"MIER_{f1}"
            sample_ls.append(new_col)

            # create new MIER column per-sample that is 1 (pass), 2 (fail), or 3 (unknown)
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

            # print(f"{f1} {df[new_col].value_counts().sort(new_col)['count'].to_list()}")
            o.write(f"{f1}\t{correct}\t{incorrect}\t{unknown}\t{(correct/(correct+incorrect)):.3f}\t{(correct/(correct+incorrect+unknown)):.3f}\n")

    if outfile:
        # write a tsv of the input vcf file coords with only the relevant samples present
        df = df.select(sample_ls)
        df = pl.concat([df_coords, df], how="horizontal")
        df.write_csv(outfile, separator="\t")


def get_parentage(parentage):
    named_f1_dt = {}

    with open(parentage) as f:
        for line in f:
            f1, p1, p2 = line.rstrip().split()
            named_f1_dt[f1] = (p1, p2)

    return named_f1_dt


if __name__ == "__main__":
    main()

