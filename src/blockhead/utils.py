#!/usr/bin/env python3

import argparse

def parse_user_input():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-p', '--parentage', type=str, required=True,
                        help='tsv parentage file')

    parser.add_argument('-v', '--vcf', type=str, required=True,
                        help='biallelic only vcf file')

    parser.add_argument('-r', '--results', type=str, required=True,
                        help='tsv summary output file')

    parser.add_argument('-o', '--outfile', type=str, required=False,
                        default=False, help='tsv file per variant')

    parser.add_argument('-t', '--threads', type=int, required=False,
                        default=1, help='max polars threads')

    parser.add_argument('-x', '--threshold', type=float, required=False,
                        default=.8, help='threshold for MIER trios to keep variant')

    parser.add_argument('-b', '--blockmode', action='store_true',
                        help='perform haplotype block functionality')

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