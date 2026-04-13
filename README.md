# BLOCKHEAD

Blockhead was designed for the analysis of trios in "Benchmarking SNP-Calling Accuracy Against Known Citrus Pedigrees Reveals Pangenome Advantages Over Linear References". It's functionality for calculated vcf-specific trio accuracy using MIER should apply to all datasets with trio data. It's advanced functionality to determine haplotype blocks `blockmode` is currently only available to third generation samples where all first (grandparental) and second generational (F1) lines are included.

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

Blockhead expects a **biallelic-SNP-only VCF**. Multi-allelic sites and indels are skipped automatically.

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

**Additional output:** one PNG per chromosome — `<chrom>_<N>_haplotypes.png`

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

**Additional output:** `<N>_haplotypes_assess_blocks.tsv`

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

## All options

```
usage: blockhead [-h] -p PARENTAGE -v VCF -d OUTDIR [-o] [-t THREADS]
                 [-n NON_MISSING] [-x THRESHOLD] [-b] [-s SMOOTH]
                 [-a ASSESS] [-w] [-c COLORS] [-k BREAKS]

options:
  -h, --help                    show this help message and exit
  -p, --parentage PARENTAGE     tsv parentage file
  -v, --vcf VCF                 biallelic only SNP input vcf file
  -d, --outdir OUTDIR           output directory for all files
  -o, --outvcf                  write a VCF of threshold-filtered variants
  -t, --threads THREADS         worker threads; requires tabix index for >1
  -n, --non_missing NON_MISSING proportion of trios with non-missing calls
                                required to keep a variant (0–1, default 0)
  -x, --threshold THRESHOLD     proportion of trios with correct calls
                                required to keep a variant (0–1, default 0)
  -b, --blockmode               perform haplotype block inference
  -s, --smooth SMOOTH           median filter half-window size in SNPs
  -a, --assess ASSESS           truth TSV for block assessment
  -w, --wrong_calls             summarize MIER-incorrect calls by type
  -c, --colors COLORS           tab-separated sample-to-color file
  -k, --breaks BREAKS           window size (bp) for wrong call binning
                                (default 1000000)
```
