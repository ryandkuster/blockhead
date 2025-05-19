#uv run python \
#    -m blockhead \
#    -v tests/data/adv_test_1.vcf \
#    -p tests/data/adv_test_1_parentage.tsv \
#    -r tests/data/adv_test_1_mier.tsv \
#    -o tests/data/adv_test_1_mier.vcf.gz \
#    -x 1 \
#    -t $1

uv run python \
    -m blockhead \
    -v tests/data/sanity_mendelian_biallelic.vcf \
    -p tests/data/mier_test_1_parentage.tsv \
    -r sanity_mendelian_biallelic_summary.tsv \
    -o sanity_mendelian_biallelic_mier.vcf \
    -x .8 \
    -t 12
