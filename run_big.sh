#uv run python \
#    -m blockhead \
#    -v /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/snp_final_scaffold_01.vcf.gz \
#    -p /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/advanced_parentage_test.tsv \
#    -r /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/adv_test_1_mier.tsv \
#    -o /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/adv_test_1_mier_rows.tsv \
#    -x .9 \
#    -t $1

for i in /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/Scaffold_0*.vcf.gz ; do
    base_name=$(basename ${i%%_snp_final.vcf.gz})
    echo $base_name
    uv run python \
        -m blockhead \
        -v $i \
        -p /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/advanced_parentage_test.tsv \
        -r /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/adv_${base_name}_MIER_summary.tsv \
        -x .9 \
        -t $1
done
