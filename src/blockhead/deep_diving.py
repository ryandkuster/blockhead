import os

import polars as pl


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