#!/usr/bin/env python3

import argparse

def parse_user_input():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-p', '--parentage', type=str, required=True,
                        help='tsv parentage file')

    parser.add_argument('-v', '--vcf', type=str, required=True,
                        help='biallelic only SNP input vcf file')

    parser.add_argument('-d', '--outdir', type=str, required=True,
                        help='output directory for all files')

    parser.add_argument('-o', '--outvcf', action='store_true',
                        help='vcf file of -x filtered variants')

    parser.add_argument('-t', '--threads', type=int, required=False,
                        default=1, help='max polars threads')

    parser.add_argument('-x', '--threshold', type=float, required=False,
                        default=.8, help='percent of children in trios with correct calls to keep variant')

    parser.add_argument('-b', '--blockmode', action='store_true',
                        help='perform haplotype block functionality')

    parser.add_argument('-c', '--colors', type=str, required=False,
                        help='optional colors file for haplotype blocks')

    args = parser.parse_args()

    return args


def gzip_test(in1):
    try:
        with open(in1) as f:
            f.readline()
        compressed = False
    except UnicodeDecodeError:
        compressed = True
    except IsADirectoryError:
        return None
    return compressed
