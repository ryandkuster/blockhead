#!/usr/bin/env python3

import argparse
import math
import numpy as np
import random

np.set_printoptions(linewidth=200)


"""
TODO list:
- [ ] add arg parsing for input
- [x] create a random number of breaks with break_no being the max
- [x] create a way to introduce random errors into the created samples
- [x] introduce window approach that uses actual base position, not just
      position in array index
- [x] create an index of positions and randomly remove segments to
      capture windows will contain varying genotypes
- [ ] add a sensitivity check for windows with few observations
- [ ] create more of a sliding window to catch more breaks in cases
      where they are close
"""


def parse_args():
    parser = argparse.ArgumentParser(description='simulate haplotype resolution')

    parser.add_argument('-s', '--snps', type=int, required=True,
                        help='number of SNPs')
    parser.add_argument('-n', '--sample_no', type=int, required=True,
                        help='number of samples to simulate')
    parser.add_argument('-b', '--break_no', type=int, required=True,
                        help='max number of breaks to simulate per sample')
    parser.add_argument('-w', '--window', type=int, required=True,
                        help='window size to use')
    parser.add_argument('-e', '--error_rate', type=float, required=True,
                        help='error rate as float')
    parser.add_argument('-c', '--spacing', type=int, required=True,
                        help='simulate spacing _ bases between SNPs')
    parser.add_argument('-p', '--pos_ls', type=str, default=None,
                        help='position list per-line text file')
    args = parser.parse_args()

    return args


def get_pos_ls(args):
    pos_ls = []
    with open(args.pos_ls) as f:
        # pos_ls = [int(line.rstrip()) for line in f]
        for line in f:
            line = line.rstrip()
            pos_ls.append(int(line))
    return pos_ls


def create_haplotypes(snps):
    """
    Create a list of random values between 0 and 2 to represent SNP calls
    """
    hap_1 = [random.randint(0, 2) for _ in range(snps)]
    hap_2 = []
    for i in hap_1:
        if i == 0:
            hap_2.append(random.choices([1, 2], k=1)[0])
        elif i == 1:
            hap_2.append(random.choices([0, 2], k=1)[0])
        elif i == 2:
            hap_2.append(random.choices([0, 1], k=1)[0])

    return [hap_1, hap_2]


def create_samples(snps, break_no, haps):
    """
    Randomly choose a starting haplotype for each sample.
    Add break points at random locations in the haplotype.
    It is possible for the number of breaks to be zero.
    """
    sample_hap = []
    start_hap = random.sample(range(0, 2), 1)[0]
    break_max = random.sample(range(0, break_no+1), 1)[0]

    # if no breaks, manually return the haplotype and empty list of breaks
    if break_max == 0:
        return haps[start_hap], []

    break_ls = sorted(random.sample(range(0, snps), break_max))
    begin_break = 0

    for end_break in break_ls:
        sample_hap += haps[start_hap][begin_break: end_break]
        begin_break = end_break
        start_hap = 0 if start_hap == 1 else 1

    if break_ls[-1] != len(haps[0]):
        sample_hap += haps[start_hap][begin_break: len(haps[0])+1]

    return sample_hap, break_ls


def introduce_error(sample_hap, error_rate):
    new_hap = []
    for gt in sample_hap:
        if random.random() <= error_rate:
            match gt:
                case 0:
                    new_hap.append(random.choices([1, 2], k=1)[0])
                case 1:
                    new_hap.append(random.choices([0, 2], k=1)[0])
                case 2:
                    new_hap.append(random.choices([0, 1], k=1)[0])
        else:
            new_hap.append(gt)

    return new_hap


def gt_to_binary(sample_array):
    """
    Make a new array (new_arrary) converting 0-2 coded types to binary.
    Use the first sample to define group membership comparisons.
    """
    bin_array = []
    for sample in sample_array:
        new_sample = []
        for idx, gt in enumerate(sample):
            if gt == sample_array[0][idx]:
                new_sample.append(0)
            else:
                new_sample.append(1)
        bin_array.append(new_sample)

    return np.array(bin_array)


def window_k_means(bin_array, window):
    """
    Iterate over the array of n samples, investigating defined windows
    """
    ref_len = len(bin_array[0])
    window_array = np.empty([math.ceil(ref_len/window), bin_array.shape[0]])

    for i in range(0, math.ceil(ref_len/window)):
        begin = window*i
        end = window*(i+1)-1
        chunk = bin_array[:, begin:end+1]
        members_ls = window_membership(chunk, window)
        window_array[i] = members_ls

    return window_array


def window_k_means_pos(bin_array, window, spacing, pos_ls):
    """
    Iterate over the array of n samples, investigating defined windows
    across uneven positions.
    The pos_idx variable simulates gt calls at random positions in the
    chromosome.
    """
    if not pos_ls:
        # pos_ls simulates missingness in gt calls based on how spaced they might be
        ref_len = len(bin_array[0])
        pos_ls = random.sample(range(1, ref_len * spacing), ref_len)
        pos_ls = sorted(pos_ls)

    max_pos = max(pos_ls)

    # list of window sizes
    window_len_ls = []

    # create empty array to fill with memberships per window
    window_array = np.empty([math.ceil(max_pos/window), bin_array.shape[0]])

    for i in range(0, math.ceil(max_pos/window)):
        # walk by chunks defined by simulated gt positions, not bin_array
        begin = window*i
        end = window*(i+1)-1
        pos_idx = [idx for idx, x in enumerate(pos_ls) if begin <= x <= end]
        chunk = bin_array[:, pos_idx]
        window_len_ls.append(len(pos_idx))
        members_ls = window_membership(chunk, len(pos_idx))
        window_array[i] = members_ls

    return window_array, window_len_ls


def window_membership(chunk, window):
    """
    Compare all samples to the first to determine membership.
    1 means same, 0 means different.
    Greater than 50% gts with ref yields ref class.
    """
    ref = chunk[0]
    members_ls = [1]

    for sample in chunk[1:, :]:
        distance = np.count_nonzero(sample != ref)
        if distance > window/2:
            members_ls.append(0)
        else:
            members_ls.append(1)

    return members_ls


def compare_membership(members_ls):
    """
    Compare membership in this window vs. last window.
    Given first sample may be at a breakpoint, this step accounts for
    sudden changes in membership.
    The reference is uninformative, so only the changes inform blocks.
    If all samples change suddenly, recode ref to 0 and swap all
    identities in remaining array.
    """
    sample_no = members_ls.shape[1] - 1
    threshold = sample_no * .5  #TODO percent of samples
    swap_ls = []
    last_row = None

    for idx, row in enumerate(members_ls):
        row = row[1:]
        if idx != 0:
            diff = np.sum(row != last_row)
            if diff >= threshold:
                swap_ls.append(idx)
        last_row = row

    return swap_ls


def swap_labels(windows_array, swap_ls):
    """
    Based on the predicted break points in swap_ls, correct the labels
    in the break ranges.
    """
    for idx in range(math.ceil(len(swap_ls)/2)):
        idx = idx*2
        try:
            windows_array[swap_ls[idx]:swap_ls[idx+1]] = 1 - windows_array[swap_ls[idx]:swap_ls[idx+1]]
        except IndexError:
            windows_array[swap_ls[idx]:] = 1 - windows_array[swap_ls[idx]:]

    return windows_array


def get_swap_pos(swap_ls, window_len_ls):
    pred_ls = []
    running_total = 0

    for idx, i in enumerate(window_len_ls):
        running_total += i
        if idx in swap_ls:
            pred_ls.append(running_total)

    return pred_ls


def main():
    args = parse_args()

    if args.pos_ls:
        pos_ls = get_pos_ls(args)
        args.snps = len(pos_ls)
    else:
        pos_ls = None

    print("simulating sample haplotypes")
    haps = create_haplotypes(args.snps)
    # haps = [[0 for i in range(0, args.snps)],[1 for i in range(0, args.snps)]] # keep for sanity of break creations

    sample_array = []
    for sample in range(0, args.sample_no):
        sample_hap, break_ls = create_samples(args.snps, args.break_no, haps)
        if sample == 0:
            truth_ls = break_ls
        sample_hap = introduce_error(sample_hap, args.error_rate)
        sample_array.append(sample_hap)

    bin_array = gt_to_binary(sample_array)
    # windows_array = window_k_means(bin_array, args.window)

    print("estimating genotypes in windows")
    windows_array, window_len_ls = window_k_means_pos(bin_array, args.window, args.spacing, pos_ls)

    print("calculating haplotypes")
    swap_ls = compare_membership(windows_array)
    windows_array = swap_labels(windows_array, swap_ls)
    pred_ls = get_swap_pos(swap_ls, window_len_ls)

    print(windows_array.T)

    print(f"predicted :   {pred_ls}")
    print(f"ground truth: {truth_ls}")


if __name__ == "__main__":
    main()
