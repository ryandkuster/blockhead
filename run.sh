uv run python \
    -m blockhead \
    -v tests/data/adv_test_1.vcf \
    -p tests/data/adv_test_1_parentage.tsv \
    -r tests/data/adv_test_1_mier.tsv \
    -o tests/data/adv_test_1_mier.vcf.gz \
    -x 1 \
    -t $1

uv run python \
    -m blockhead \
    -v tests/data/mier_test_1.vcf \
    -p tests/data/mier_test_1_parentage.tsv  \
    -r summary.delete.txt \
    -o delete_test.vcf.gz \
    -x .7 \
    -t $1
