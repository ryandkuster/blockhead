OUTDIR=output/adv_assess
mkdir -p $OUTDIR

blockhead \
    --parentage input/adv_test_1_parentage.tsv \
    --vcf input/adv_test_1.vcf.gz \
    --outdir $OUTDIR \
    --threads 12 \
    --non_missing 0 \
    --threshold 0 \
    --blockmode \
    --assess output/adv_smooth/2_haplotypes_5000_smooth_blocks.tsv \
    --colors input/adv_test_1_colors.tsv

