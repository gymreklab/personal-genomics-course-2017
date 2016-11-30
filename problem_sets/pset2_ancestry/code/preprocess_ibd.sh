#!/bin/bash

INVCF=/oasis/projects/nsf/csd524/mgymrek/data/ps2/ALL.chr16.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz
SUBVCF=/oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_ibd.vcf.gz
PREFIX=/oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_ibd

# Get subset of samples
zcat ${INVCF} | \
    vcf-subset -c NA19713,NA19985,NA19625,NA20414,NA19660 | \
    bgzip -c > ${SUBVCF}
tabix -p vcf ${SUBVCF}

# Conver to Beagle format
zcat ${SUBVCF} | \
    java -jar -Xmx512m /home/mgymrek/sources/vcf2beagle.jar -9 ${PREFIX}

