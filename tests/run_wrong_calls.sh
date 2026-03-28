OUTDIR=output/adv_wrong_calls
mkdir -p $OUTDIR

blockhead \
    --parentage input/adv_test_1_parentage.tsv \
    --vcf input/adv_test_1.vcf.gz \
    --outdir $OUTDIR \
    --threads 12 \
    --non_missing 0 \
    --threshold 0 \
    --wrong_calls \
    --breaks 10_000

