## Installation

Install UV for an easy experience:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Run the following from blockhead repo to make a package available on command line.

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

```bash
uv run python \
    -m blockhead \
    -v tests/data/adv_test_1.vcf \
    -p tests/data/adv_test_1_parentage.tsv \
    -r tests/data/adv_test_1_mier.tsv \
    -o tests/data/adv_test_1_mier.vcf.gz \
    -x 1 \
    -t 12
```

```
options:
  -h, --help            show this help message and exit
  -p PARENTAGE, --parentage PARENTAGE
                        tsv parentage file
  -v VCF, --vcf VCF     biallelic only vcf file
  -r RESULTS, --results RESULTS
                        tsv summary output file
  -o OUTFILE, --outfile OUTFILE
                        tsv file per variant
  -t THREADS, --threads THREADS
                        max polars threads
  -x THRESHOLD, --threshold THRESHOLD
                        threshold for MIER trios to keep variant
  -b, --blockmode       perform haplotype block functionality
```

