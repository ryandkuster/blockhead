## Installation

Install UV for an easy experience:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Run the following while in the blockhead repo to make a package available on command line.

```bash
uv venv --python 3.12
```

...then:
```bash
pip install -e .
```

...if on macOS, try:
```bash
pipx install -e .
```

## Usage

```bash
blockhead \
    -v tests/data/adv_test_1.vcf \
    -p tests/data/adv_test_1_parentage.tsv \
    -r tests/data/adv_test_1_mier.tsv \
    -o tests/data/adv_test_1_mier.vcf.gz \
    -x 1 \
    -t 12
```

For more help, use:
```bash
blockhead --help
```

```
usage: blockhead [-h] -p PARENTAGE -v VCF -d OUTDIR [-o] [-t THREADS] [-x THRESHOLD] [-b]
                 [-c COLORS]

options:
  -h, --help            show this help message and exit
  -p, --parentage PARENTAGE
                        tsv parentage file
  -v, --vcf VCF         biallelic only SNP input vcf file
  -d, --outdir OUTDIR   output directory for all files
  -o, --outvcf          vcf file of -x filtered variants
  -t, --threads THREADS
                        max polars threads
  -x, --threshold THRESHOLD
                        percent of children in trios with correct calls to keep variant
  -b, --blockmode       perform haplotype block functionality
  -c, --colors COLORS   optional colors file for haplotype blocks
```

