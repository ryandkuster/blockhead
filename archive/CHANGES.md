# Changes vs. Original Implementation

The original code lives in `archive/blockhead/`. The key driver for the rewrite
was memory: the original loaded the entire VCF into a Polars DataFrame at
startup, which became impractical for large files.

---

## 1. VCF loading strategy

**Original:** `pl.read_csv()` on the entire VCF at startup. The full file sat
in memory as a string-typed Polars DataFrame for the duration of the run.

**New:** `cyvcf2` streaming via `stream_chromosomes()` (sequential) or
`fetch_chrom()` (parallel, requires tabix index). One chromosome's worth of
data is held in memory at a time and released before the next is loaded.

---

## 2. Genotype encoding

**Original:** Genotypes were recoded as strings `"0"`, `"1"`, `"2"` via regex
substitution in `recode_vcf()`. All subsequent filtering compared string
columns.

**New:** `cyvcf2` exposes `variant.gt_types` directly as a uint8 array using
its own encoding: `HOM_REF=0`, `HET=1`, `UNKNOWN=2`, `HOM_ALT=3`. Genotypes
are stored in a numpy matrix (`gts_matrix`) of shape
`(n_variants, n_samples)`, dtype uint8. No string conversion anywhere in the
hot path.

---

## 3. MIER computation

**Original:** Computed per-F1 MIER as a new string column on the full
DataFrame using chained Polars `when/then` expressions comparing string
genotypes.

**New:** `compute_mier()` in `df_manip.py` operates on three uint8 numpy
column slices and returns an int8 array. Called once per F1 per chromosome;
results stored in `mier_dict` keyed by F1 name.

---

## 4. Per-sample haplotype inference

**Original:** `infer_haplotypes()` received the full per-chromosome Polars
DataFrame (already in memory) and applied four sequential `.filter()` passes
plus two `.with_columns()` to produce a result DataFrame. One pre-built
DataFrame was reused across all advanced hybrids on a chromosome.

**New:** `infer_haplotypes()` works entirely on numpy arrays. It receives
`gts_matrix`, `mier_dict`, and a `lin_dt` dict describing the lineage, pulls
the six relevant column slices directly, builds a single combined boolean mask
in one pass, and returns `(pos_filtered, parent_origin)` as numpy arrays. No
Polars DataFrame is constructed per sample. The `x` and `y` lists fed into the
scatter plot and smoothing code are identical in content to the original.

---

## 5. Parallelism

**Original:** Fully sequential — one chromosome after another in a single
pass.

**New:** If the VCF has a tabix index, chromosomes are dispatched to a
`ProcessPoolExecutor` (one worker per chromosome, up to `--threads`). Each
worker calls `fetch_chrom()` independently. Without a tabix index the code
falls back to sequential streaming with a warning.

---

## 6. MIER filtering threshold

**Original:** Filtering was applied to the full in-memory DataFrame using
Polars `percent_mier_correct` and `percent_non_missing` columns.

**New:** Per-chromosome, a numpy stack of all MIER arrays is built and the
pass mask is computed vectorially. Passing `(chrom, pos)` pairs are collected
into a set and used in a second streaming pass to write the filtered VCF if
`--outvcf` is requested.

---

## What did not change

- The `median_filter` smoothing logic (`--smooth`) is unchanged from the
  original pure-Python implementation.
- The `assess_blocks` logic (`--assess`) is unchanged.
- All matplotlib figure layout, axes structure, scatter parameters, and output
  filenames are identical — output images are pixel-equivalent to the original.
- The parentage file format and `parentage.py` logic are unchanged.
- The `--wrong_calls` analysis path is unchanged in structure.
