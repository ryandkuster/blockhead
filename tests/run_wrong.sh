#blockhead \
#    -v data/sanity_mendelian_biallelic.vcf \
#    -p data/mier_test_1_parentage.tsv \
#    -d output/mier_wrong \
#    -o \
#    -w \
#    -x .8 \
#    -t 12

blockhead \
    -v data/adv_test_1.vcf.gz \
    -p data/adv_test_1_parentage.tsv \
    -d output/adv_wrong \
    -c data/adv_test_2_wrong_colors.tsv \
    -o \
    -x .5 \
    -w \
    -k 100_000 \
    -t 12

