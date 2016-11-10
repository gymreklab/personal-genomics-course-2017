#!/bin/bash

# TODO do for whole file

VCFFILE=/Users/gymrek/Downloads/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

zcat < $VCFFILE | vcf-to-tab | head -n 100000 | ./get_vcf_matrix.py ../data/ps2_pca

