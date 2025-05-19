uv run python \
    -m blockhead \
    -v /Users/ryankuster/Downloads/tmp_citrus_parentage/snp_final_bcftools.vcf.gz \
    -p /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/citrus_parentage_bcftools.tsv \
    -r /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/summaries/bcftools_40_mier.tsv \
    -o /Users/ryankuster/Downloads/tmp_citrus_advanced_parentage/vcf_mier/bcftools_40_mier.vcf \
    -x .8 \
    -t 12
