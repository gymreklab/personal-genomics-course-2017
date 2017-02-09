#!/bin/bash

NIST=/oasis/projects/nsf/csd524/mgymrek/data/ps4/NIST_NA12878_hg38_chr22.tab
SNPCALLS=/oasis/projects/nsf/csd524/$USER/ps4/results/NA12878_chr22_snps.tab

# Only look at non homref - otherwise this file is huge
cat ${SNPCALLS} | awk '($4!= ".")' > ${SNPCALLS}.onlysnps

# Load to compare script
./ps4_comparesnps.py \
    --nist ${NIST} \
    --snpcalls ${SNPCALLS}.onlysnps


