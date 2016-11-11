#!/bin/bash

VCFFILE=ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# Extract only SNPs used in 23andme
intersectBed -header -a ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -b 23andme_snps_hg19.bed | \
    bgzip -c > ps2_pca.genotypes.vcf.gz
tabix -p vcf ps2_pca.genotypes.vcf.gz
