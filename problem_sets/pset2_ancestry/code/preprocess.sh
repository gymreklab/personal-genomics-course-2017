#!/bin/bash

VCFFILE=ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# Extract only SNPs used in 23andme
zcat ${VCFFILE} | grep -v "VT=SV" | \
    intersectBed -header -a stdin -b 23andme_snps_hg19.bed.gz | \
    bgzip -c > ps2_pca.genotypes.vcf.gz
tabix -p vcf ps2_pca.genotypes.vcf.gz

# Get matrix
zcat ps2_pca.genotypes.vcf.gz | vcf-to-tab | get_vcf_matrix.py ps2_pca