# BLOCKHEAD

<img src=assets/images/blockhead.png align="left" width="250" alt="Blockhead Logo">

Blockhead was designed for the analysis of trios in "Benchmarking SNP-Calling Accuracy Against Known Citrus Pedigrees Reveals Pangenome Advantages Over Linear References". It's functionality for calculated vcf-specific trio accuracy using MIER should apply to all datasets with trio data. It's advanced functionality to determine haplotype blocks `blockmode` is currently only available to third generation samples where all first (grandparental) and second generational (F1) lines are included.

If you have any questions on the use of or desired functionality for blockhead, please contact Ryan Kuster (rkuster@utk.edu).

## Installation

Install UV for an easy experience:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Run the following while in the blockhead repo to create and activate a virtual environment:

```bash
uv venv --python 3.12
source .venv/bin/activate
```

Then install blockhead as an editable package:

```bash
pip install -e .
```

On macOS, `pipx` is an alternative that manages the environment automatically:

```bash
pipx install -e .
```

## Input files

### VCF

Blockhead expects a VCF with at most one ALT allele per site. Multi-allelic sites (comma in the ALT field) are skipped automatically; indels and longer sequences are retained.

For parallel execution (`--threads > 1`) the VCF must have an index (`.csi`). Without one, blockhead falls back to single-threaded streaming with a warning. If compressed, only bgzipped inputs are handled, not gzipped.

### Parentage file (`--parentage`)

A tab-separated file with no header. Each row defines one trio:

```
<child>    <parent1>    <parent2>
```

Sample names must match the column headers in the VCF exactly. For advanced hybrid analysis (blockmode), an F1 sample can appear as a parent in a second row:

```
a_b     a_a     b_b
ab_c1   c_c     a_b
ab_c2   c_c     a_b
```

Here `a_b` is both an F1 in row 1 and a parent in rows 2–3, making `ab_c1` and `ab_c2` advanced hybrids.

#### F1 sibling crosses

Two F1s that share the same original parents can themselves be crossed. Blockhead detects this automatically when both parents of a child are F1s with identical grandparents:

```
a_b1    a_a     b_b
a_b2    a_a     b_b
a_b12   a_b1    a_b2
```

Here `a_b1` and `a_b2` are both F1s from the same `a_a × b_b` cross. Their offspring `a_b12` is an F1 sibling cross. Blockhead treats it separately from a standard advanced hybrid: informative sites require both F1 parents to be heterozygous, and child genotype is assigned to the grandparent (p1 or p2) whose homozygous state it matches. Heterozygous children are excluded as ambiguous. F1 sibling crosses produce their own set of output files (see `--blockmode` outputs below).

### Colors file (`--colors`)

A tab-separated file with no header, mapping parental (not F1) sample IDs to a matplotlib color string:

```
a_a    "orange"
b_b    "green"
a_b    "purple"
```

Required when `--blockmode` is used.

## Usage

### Basic MIER summary

Compute Mendelian Inheritance Error Rate (MIER) for all trios and write summary tables:

```bash
blockhead \
    --parentage tests/input/adv_test_1_parentage.tsv \
    --vcf tests/input/adv_test_1.vcf.gz \
    --outdir output/ \
    --threshold 0.9 \
    --non_missing 0.8
```

**Outputs:**
- `MIER_summary.tsv` — per-F1 correct/incorrect/unknown counts and ratios
- `MIER_summary_homozygous_parents.tsv` — same, restricted to sites where parents are homozygous-opposite

### Filtered VCF (`--outvcf`)

Add `--outvcf` to write a VCF containing only variants that pass the `--threshold` and `--non_missing` filters:

```bash
blockhead \
    --parentage tests/input/adv_test_1_parentage.tsv \
    --vcf tests/input/adv_test_1.vcf.gz \
    --outdir output/ \
    --threshold 0.9 \
    --non_missing 0.8 \
    --outvcf
```

**Additional output:** `MIER_filtered.vcf`

### Haplotype block plots (`--blockmode`)

Infer and plot parental haplotype origins for advanced hybrids. Requires a `--colors` file:

```bash
blockhead \
    --parentage tests/input/adv_test_1_parentage.tsv \
    --vcf tests/input/adv_test_1.vcf.gz \
    --outdir output/ \
    --threads 12 \
    --threshold 0 \
    --non_missing 0 \
    --blockmode \
    --colors tests/input/adv_test_1_colors.tsv
```

**Additional outputs:**
- `<chrom>_<N>_haplotypes.png` — haplotype plot per chromosome (standard advanced hybrids)
- `<chrom>_<N>_f1_cross_haplotypes.png` — haplotype plot per chromosome (F1 sibling crosses, when present)

### Smoothing (`--smooth`)

Apply a median filter over haplotype block assignments. `--smooth` takes an integer that sets the half-window size in number of SNPs:

```bash
blockhead \
    --parentage tests/input/adv_test_1_parentage.tsv \
    --vcf tests/input/adv_test_1.vcf.gz \
    --outdir output/ \
    --threads 12 \
    --blockmode \
    --smooth 5000 \
    --colors tests/input/adv_test_1_colors.tsv
```

**Additional outputs:**
- `<chrom>_<N>_haplotypes_smooth.png` — smoothed haplotype plot per chromosome
- `<chrom>_<N>_haplotypes_breaks.png` — breakpoint histogram per chromosome
- `<N>_haplotypes_<smooth>_smooth_blocks.tsv` — per-position smoothed block assignments across all chromosomes
- `<chrom>_<N>_f1_cross_haplotypes_smooth.png` — smoothed haplotype plot for F1 sibling crosses (when present)
- `<chrom>_<N>_f1_cross_haplotypes_breaks.png` — breakpoint histogram for F1 sibling crosses (when present)
- `<N>_f1_cross_haplotypes_<smooth>_smooth_blocks.tsv` — smoothed block assignments for F1 sibling crosses (when present)

### Block assessment (`--assess`)

Compare smoothed block calls against a truth set. The truth file has the same format as the smooth blocks TSV output (typically the output of a previous `--smooth` run used as ground truth):

```bash
blockhead \
    --parentage tests/input/adv_test_1_parentage.tsv \
    --vcf tests/input/adv_test_1.vcf.gz \
    --outdir output/ \
    --threads 12 \
    --blockmode \
    --assess output/adv_smooth/2_haplotypes_5000_smooth_blocks.tsv \
    --colors tests/input/adv_test_1_colors.tsv
```

**Additional output:**
- `<N>_haplotypes_assess_blocks.tsv` — block accuracy assessment for standard advanced hybrids
- `<N>_f1_cross_haplotypes_assess_blocks.tsv` — block accuracy assessment for F1 sibling crosses (when present)

### Wrong call analysis (`--wrong_calls`)

Summarize MIER-incorrect calls by type and genomic window. This mode exits after writing its output and skips all other analysis:

```bash
blockhead \
    --parentage tests/input/adv_test_1_parentage.tsv \
    --vcf tests/input/adv_test_1.vcf.gz \
    --outdir output/ \
    --threads 12 \
    --wrong_calls \
    --breaks 10000
```

**Outputs:**
- `MIER_per_snp_summary.tsv` — per-SNP MIER values for all F1s
- `<windows>_windows_wrong_calls.tsv` — wrong call counts binned by genomic window
- `<chrom>_<N>_<window>_wrong_calls.png` — wrong call plot per chromosome

### Quality filtering (`--quality`)

Skip variants whose QUAL score falls below a minimum value. By default no quality filter is applied:

```bash
blockhead \
    --parentage tests/input/adv_test_1_parentage.tsv \
    --vcf tests/input/adv_test_1.vcf.gz \
    --outdir output/ \
    --quality 30
```

## All options

```
usage: blockhead [-h] -p PARENTAGE -v VCF -d OUTDIR [-o] [-t THREADS]
                 [-n NON_MISSING] [-x THRESHOLD] [-q QUALITY] [-b] [-s SMOOTH]
                 [-a ASSESS] [-w] [-c COLORS] [-k BREAKS]

options:
  -h, --help                    show this help message and exit
  -p, --parentage PARENTAGE     tsv parentage file
  -v, --vcf VCF                 input vcf file
  -d, --outdir OUTDIR           output directory for all files
  -o, --outvcf                  write a VCF of threshold-filtered variants
  -t, --threads THREADS         worker threads; requires tabix index for >1
  -n, --non_missing NON_MISSING proportion of trios with non-missing calls
                                required to keep a variant (0–1, default 0)
  -x, --threshold THRESHOLD     proportion of trios with correct calls
                                required to keep a variant (0–1, default 0)
  -q, --quality QUALITY         minimum QUAL score; variants below this value
                                are skipped (default: no filter)
  -b, --blockmode               perform haplotype block inference
  -s, --smooth SMOOTH           median filter half-window size in SNPs
  -a, --assess ASSESS           truth TSV for block assessment
  -w, --wrong_calls             summarize MIER-incorrect calls by type
  -c, --colors COLORS           tab-separated sample-to-color file
  -k, --breaks BREAKS           window size (bp) for wrong call binning
                                (default 1000000)
```
