blockhead \
    -v data/sanity_mendelian_biallelic.vcf \
    -p data/mier_test_1_parentage.tsv \
    -d output/mier \
    -o \
    -x .8 \
    -t 12

blockhead \
    -v data/adv_test_1.vcf.gz \
    -p data/adv_test_1_parentage.tsv \
    -c data/adv_test_1_colors.tsv \
    -d output/adv \
    -o \
    -x .5 \
    -b \
    -t 12
