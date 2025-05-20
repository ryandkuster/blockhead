blockhead \
    -v data/sanity_mendelian_biallelic.vcf \
    -p data/mier_test_1_parentage.tsv \
    -r output/sanity_mendelian_biallelic_summary.tsv \
    -o output/sanity_mendelian_biallelic_mier.vcf \
    -x .8 \
    -t 12

blockhead \
    -v data/mier_test_1.vcf \
    -p data/mier_test_1_parentage.tsv \
    -r output/mier_test_1_summary.tsv \
    -o output/mier_test_1_mier.vcf \
    -x .8 \
    -t 12
