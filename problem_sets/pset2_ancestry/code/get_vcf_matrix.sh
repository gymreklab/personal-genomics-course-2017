#!/bin/bash

# Replace with 23andme VCF
VCFFILE=../data/ps2_pca.genotypes.vcf.gz
zcat $VCFFILE | vcf-to-tab | ./get_vcf_matrix.py ../data/ps2_pca

