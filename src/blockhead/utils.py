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

    parser.add_argument('-n', '--non_missing', type=float, required=False,
                        default=0.0, help='proportion of children in trios with non-missing calls to keep variant (max 1)')

    parser.add_argument('-x', '--threshold', type=float, required=False,
                        default=0.0, help='proportion of children in trios with correct calls to keep variant (max 1)')

    parser.add_argument('-b', '--blockmode', action='store_true',
                        help='perform haplotype block functionality')

    parser.add_argument('-s', '--smooth', type=int, required=False,
                        default=0, help='smooth haplotype blocks using median filter')

    parser.add_argument('-a', '--assess', type=str, required=False,
                        help='use smooth haplotype blocks as ground truth vs input')

    parser.add_argument('-w', '--wrong_calls', action='store_true',
                        help='assess wrong calls only')

    parser.add_argument('-c', '--colors', type=str, required=False,
                        help='optional colors file for haplotype blocks')

    parser.add_argument('-q', '--quality', type=int, required=False,
                        default=None,
                        help='minimum variant QUAL score to retain (VCF field 6); '
                             'variants below this value are excluded')

    parser.add_argument('-k', '--breaks', type=int, required=False,
                        default=1000000, help='break size for wrong call windows')

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
