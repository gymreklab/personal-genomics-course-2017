#!/bin/bash

VCFFILE1=/oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_pca.genotypes.vcf.gz
VCFFILE2=/oasis/projects/nsf/csd524/mgymrek/23andMe/genome_Melissa_Gymrek_v3_Full_20161129092229_nochr.vcf.gz
MYFILE=/oasis/projects/nsf/csd524/mgymrek/23andMe/genome_Melissa_Gymrek_v3_Full_20161129092229.vcf.gz

zcat ${MYFILE} | sed 's/^chr//' | bgzip -c > ${VCFFILE2}
tabix -p vcf ${VCFFILE2}

vcf-merge ${VCFFILE1} ${VCFFILE2} | bgzip -c > ps2_pca.genotypes_23andme.vcf.gz
tabix -p vcf ps2_pca.genotypes_23andme.vcf.gz
zcat ps2_pca.genotypes_23andme.vcf.gz | vcf-to-tab | get_vcf_matrix.py ps2_pca_23andme

