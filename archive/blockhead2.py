#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "hmmlearn",
#     "matplotlib",
#     "polars==0.20.16",
#     "scipy",
# ]
# ///

"""
Notes:
    - Using the df object that is (optionally) output from
      mendelian_error.py, 

Inputs:
    - correct = df.filter(pl.col(new_col) == 1).height
    - incorrect = df.filter(pl.col(new_col) == 2).height
    - unknown = df.filter(pl.col(new_col) == 3).height
    - file showing the expected combos of parents to produce F1s
      (3 columns) (same used for mendelian_error script)
        - example line format: f1\tp1\tp2

Confirm f1 genotype is possible given parental genotypes:
    - 0/0 cannot originate from p1 = 1/1 OR p2 = 1/1
    - 0/1 cannot originate from (p1 = 1/1 and p2 = 1/1) OR (p1 = 0/0 and p2 = 0/0)
    - 1/1 cannot originate from p1 = 0/0 OR p2 = 0/0
    - order of 0/1 vs. 1/0 is not treated as phased
    - multiallelic genotypes (alleles = 2+) or missing (./1 etc.) in f1, p1, ori
      p2 counted as unknown
"""

import math
import matplotlib.pyplot as plt
import numpy as np
import os
import polars as pl
import sys

from itertools import combinations
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import cophenet

np.set_printoptions(linewidth=100000)
np.set_printoptions(threshold=np.inf)


def get_parent_dt(parentage):
    parent_dt = {}
    cross_ls = []

    with open(parentage) as f:
        for line in f:
            f1, p1, p2 = line.rstrip().split()
            if p1 not in parent_dt:
                parent_dt[p1] = []
            if p2 not in parent_dt:
                parent_dt[p2] = []
            parent_dt[p1].append(f1)
            parent_dt[p2].append(f1)
            cross_ls.append((p1, p2))
    
    cross_ls = list(set(cross_ls))

    return parent_dt, cross_ls


def haplotype_per_cross(parent_dt, mier_p1, mier_p2, mier_df, window, cutoff):
    """
    Designating one parent as the parent of interest (mier_p1), iterate
    one scaffold at a time to see the calls that likely were contributed
    by this parent.
    """

    chrom_ls = sorted(mier_df["#CHROM"].unique().to_list())
    for chrom in chrom_ls:

        # filter to a single chromosome
        print(f"subsetting {chrom}")
        tmp_df = mier_df.filter(
            (mier_df["#CHROM"] == chrom)
        )

        # let's just start by looking at parent 1 (p1)
        f1_ls = [f1 for f1 in parent_dt[mier_p1] if f1 in parent_dt[mier_p2]]
        f1_mier_ls = [f"MIER_{f1}" for f1 in f1_ls]

        # get the sites where p2 is homozygous (0 or 2)
        tmp_df = tmp_df.filter(
            ((tmp_df[mier_p2] == "0") | (tmp_df[mier_p2] == "2")) &
            (tmp_df[mier_p1] == "1")
        )

        # filter sites where progeny from this cross are mier correct (1)
        tmp_df = tmp_df.filter(pl.all_horizontal([tmp_df[col] == 1 for col in f1_mier_ls]))

        # filter for non-missing
        # tmp_df = tmp_df.filter(pl.all_horizontal([tmp_df[col].is_in([0, 1, 2]) for col in f1_ls]))


        # select only the necessary columns
        tmp_df = tmp_df.select(["#CHROM", "POS"] + f1_ls + [mier_p1, mier_p2])

        pos_ls = tmp_df["POS"].to_list()
        tmp_df = tmp_df.select(tmp_df.columns[2:])
        tmp_df = tmp_df.drop([mier_p1, mier_p2])

        np_haps = tmp_df.to_numpy().T

        # Remove completely uninformative sites.
        matches = []
        for i in range(np_haps.shape[1]):
            if np.unique(np_haps[:, i]).size == 1:
                matches.append(i)
        np_haps = np.delete(np_haps, matches, axis=1)
        pos_ls = [i for idx, i in enumerate(pos_ls) if idx not in matches]

        # Find positions that define haplotype regions.
        useful_positions = window_assessment(np_haps, window, pos_ls, cutoff)
        useful_positions.sort()
        print(len(useful_positions))

        # Convert the positions to indices to further filter np_haps.
        idx_ls = [idx for idx, i in enumerate(pos_ls) if i in useful_positions]
        np_haps = np_haps[:, idx_ls]

        # Compare to the first sample to get binary labels.
        bin_array = gt_to_binary(np_haps)

        # Here be dragons.
        windows_array, distance_array, window_len_ls = window_k_means_pos(bin_array, window, useful_positions)

        # Compare membership changes signaling that ref sample swapped.
        swap_ls = compare_membership(windows_array)
        windows_array = swap_labels(windows_array, swap_ls)
        pred_ls = get_swap_pos(swap_ls, window_len_ls, useful_positions)

        print(f"predicted swap points windows for ref sample : {pred_ls}")

        # Make a plot of the original distance and predictions.
        fig, axs = plt.subplots(4, 1, figsize=(30, 7), sharex=False)
        axs[0].plot(pos_ls, np.zeros_like(pos_ls), '|', markersize=10)
        axs[1].plot(useful_positions, np.zeros_like(useful_positions), '|', markersize=10)
        axs[2].imshow(distance_array.T, cmap='viridis', interpolation='nearest', aspect="auto")
        axs[3].imshow(windows_array.T, cmap='viridis', interpolation='nearest', aspect="auto")
        plt.show()


def window_assessment(np_haps, window, pos_ls, cutoff):
    """
    Iterate over the array of n samples, investigating defined windows
    across uneven positions.
    """
    useful_positions = []
    max_pos = max(pos_ls)

    # Get all pairwise combinations of list indices
    list_pairs = list(combinations(range(np_haps.shape[0]), 2))

    remove_indices = []

    # for begin in range(0, max_pos-window):
    for begin in range(0, max_pos-window, window//10):
        end = begin+window
        pos_idx = [idx for idx, x in enumerate(pos_ls) if begin <= x <= end]

        if len(pos_idx) < 10:
            remove_indices.append(begin)
            print(f"{begin} {end} : {len(pos_idx)} variants, SKIP")
            continue

        chunk = np_haps[:, pos_idx]
        chunk_pos = [x for idx, x in enumerate(pos_ls) if begin <= x <= end]

        print(f"{begin} {end} : {len(pos_idx)} variants")

        # Dictionary to store similarity scores for each pair
        similarity_data = {}

        for i, j in list_pairs:
            # Find positions where lists have the same values
            matching_positions = np.where(chunk[i] == chunk[j])[0]
    
            # Calculate similarity as the proportion of matching positions
            similarity = len(matching_positions) / len(chunk[i])
    
            # Store both the similarity score and the matching positions
            similarity_data[(i, j)] = {
                'score': similarity,
                'positions': matching_positions
            }

        # Create a distance matrix for the given chunk.
        distance_matrix = np.zeros((np_haps.shape[0], np_haps.shape[0]))

        # Based on similarity of the chunk, fill out the matrix.
        for (i, j), similarity in similarity_data.items():
            distance = 1 - similarity["score"]  # distance = 1 - similarity
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance  # symmetric

        # Cluster into groups based on t threshold with fcluster.
        condensed_distance = squareform(distance_matrix)
        Z = linkage(condensed_distance, method='average')
        clusters_1 = fcluster(Z, t=0.5, criterion='distance')  # t can be tuned
        num_clusters = len(set(clusters_1))

        # Use cophenetic distance to guage how informative the window is.
        c, coph_dists = cophenet(Z, condensed_distance)

        print(c)

        if num_clusters == 2 and c >= cutoff:
            group_1_useful = get_informative_sites(similarity_data, clusters_1, 1)
            group_2_useful = get_informative_sites(similarity_data, clusters_1, 2)
            window_useful = list(group_1_useful.intersection(group_2_useful))
            chunk_useful = [chunk_pos[i] for i in window_useful]
            print(len(chunk_useful))
            useful_positions = list(set(useful_positions + chunk_useful))
            clusters_1 = correct_blocks(clusters_1)
        else:
            remove_indices.append(begin)

    return useful_positions


def get_informative_sites(similarity_data, clusters_1, group):
    """
    Now that we have a clearly defined window with only 2 groups, find
    the positions that are commonly useful for defining both groups.
    """
    group_pos_ls = []
    group_ls = [idx for idx, i in enumerate(clusters_1) if i == group]
    list_pairs = list(combinations(group_ls, 2))

    for pair in list_pairs:
        group_pos_ls.append(similarity_data[pair]["positions"].tolist())
    return find_common_elements(group_pos_ls)


def find_common_elements(group_pos_ls):
    try:
        common_elements = set(group_pos_ls[0])
    except IndexError:
        return set([])
    for i in group_pos_ls[1:]:
        common_elements = common_elements.intersection(set(i))
    
    return set(common_elements)


def correct_blocks(clusters_1):
    if clusters_1[0] == 1:
        return clusters_1
    else:
        return np.array([1 if i == 2 else 2 for i in clusters_1])


def gt_to_binary(sample_array):
    """
    Make a new array (bin_array) converting 0-2 coded types to binary.
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


def window_k_means_pos(bin_array, window, pos_ls):
    """
    Iterate over the array of n samples, investigating defined windows
    across uneven positions.
    The pos_idx variable simulates gt calls at random positions in the
    chromosome.
    """
    # pos_ls simulates missingness in gt calls based on how spaced they might be
    ref_len = len(bin_array[0])
    max_pos = max(pos_ls)

    # list of window sizes
    window_len_ls = []

    # create empty arrays to fill with memberships per window
    window_array = np.empty([math.ceil(max_pos/window), bin_array.shape[0]])
    distance_array = np.empty([math.ceil(max_pos/window), bin_array.shape[0]])

    remove_indices = []

    for i in range(0, math.ceil(max_pos/window)):
        begin = window*i
        end = window*(i+1)-1
        pos_idx = [idx for idx, x in enumerate(pos_ls) if begin <= x <= end]

        #TODO define a threshold for missingness or muddy distance to skip.
        if len(pos_idx) == 0:
            remove_indices.append(i)
            continue

        chunk = bin_array[:, pos_idx]
        window_len_ls.append(len(pos_idx))
        members_ls, distance_ls = window_membership(chunk, len(pos_idx))
        window_array[i] = members_ls
        distance_array[i] = distance_ls

    #TODO correct the array to remove windows with low variant counts
    window_array = np.delete(window_array, remove_indices, axis=0)
    return window_array, distance_array, window_len_ls


def window_membership(chunk, window):
    """
    Compare all samples to the first to determine membership.
    0 means same, 1 means different.
    Greater than 50% gts with ref yields ref class.
    """
    ref = chunk[0]
    members_ls = [1]
    distance_ls = [0]

    for sample in chunk[1:, :]:
        distance = np.count_nonzero(sample != ref)
        distance_ls.append(distance/window)

        # if distance > window/2:
        if distance > window/3:
            members_ls.append(0)
        else:
            members_ls.append(1)

    print(distance_ls)
    return members_ls, distance_ls


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


def plot_haplotypes(result, f1_ls, p1, p2, chrom, out_dir):
    f1_len = len(f1_ls)

    fig, axes = plt.subplots(nrows=f1_len, sharex=True, figsize=(30, f1_len/2), gridspec_kw={'hspace': 0.3})

    for idx, f1 in enumerate(f1_ls):
        label_name = "_".join(f1.split("_")[1:4])

        x = result["POS"].to_list()
        y = result[f"PHASE_{f1}"].to_list()

        colors_ls = ["#008080" if i == 1 else "#FFA07A" for i in y]
        axes[idx].vlines(x, ymin=0, ymax=.05, colors=colors_ls, linewidth=0.1, alpha=1)
        axes[idx].set_ylabel(f"{label_name}", labelpad=40, loc="center", rotation=0)
        axes[idx].set_yticks([])

        axes[-1].set_xlabel("Position (10Mb)")  # Only label x-axis on bottom plot
        fig.suptitle(f"{chrom} {p1}", fontsize=12, y=0.95)

    plt.savefig(os.path.join(out_dir, f"{chrom}_{p1}_x_{p2}_parental_haplotypes.png"))


def get_swap_pos(swap_ls, window_len_ls, pos_ls):
    """
    Get positions in bp of where breaks (swaps) are occurring.
    Each window contains a known number of entries originating from the 
    original vcf. We grab the associated positions from the vcf in pos_ls,
    so indexing at the swap points will translate into original vcf terms.
    """
    pred_ls = []
    running_total = 0

    # print(window_len_ls) #TODO
    # print(swap_ls) #TODO

    for idx, i in enumerate(window_len_ls):
        running_total += i
        # print(running_total) #TODO
        if idx in swap_ls:
            pred_ls.append(pos_ls[running_total-1])

    return pred_ls


def main():
    parentage = sys.argv[1]     # parentage file
    mier = sys.argv[2]          # mendelian_error output df (tsv of gt/mier)
    window = int(sys.argv[3])   # window size
    cutoff = float(sys.argv[4]) # cophenetic similarity minimum

    parent_dt, cross_ls = get_parent_dt(parentage)

    print("opening mier file")
    mier_df = pl.read_csv(mier,
                     separator="\t",
                     null_values="..")

    for mier_p1, mier_p2 in cross_ls:
        print(f"parent 1 : {mier_p1}")
        print(f"parent 2 : {mier_p2}")
        haplotype_per_cross(parent_dt, mier_p1, mier_p2, mier_df, window, cutoff)

    for mier_p2, mier_p1 in cross_ls:
        print(f"parent 1 : {mier_p1}")
        print(f"parent 2 : {mier_p2}")
        haplotype_per_cross(parent_dt, mier_p1, mier_p2, mier_df, window, cutoff)


if __name__ == "__main__":
    main()
