uv run python \
    -m blockhead \
    -v /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/citrus_bcftools_combined_SNP_PASS_Q30_BIALLELIC.vcf.gz \
    -p /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/citrus_parentage_bcftools.tsv \
    -r /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/summaries/bcftools_mier.tsv \
    -o /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/vcf_mier/bcftools_mier.vcf.gz \
    -x .8 \
    -t 12
