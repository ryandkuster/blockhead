OUTDIR=output/adv_blockmode
mkdir -p $OUTDIR

blockhead \
    --parentage input/adv_test_1_parentage.tsv \
    --vcf input/adv_test_1.vcf.gz \
    --outdir $OUTDIR \
    --threads 12 \
    --non_missing 0 \
    --threshold 0 \
    --blockmode \
    --colors input/adv_test_1_colors.tsv

