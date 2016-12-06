B1;2c#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 1000
#SBATCH --get-user-env

DATADIR=/oasis/projects/nsf/csd524/mgymrek/data/ps2/
PREFIX=${DATADIR}/ps2_impute
SAMPLE=NA20340

# Get all SNPs for one LWK sample
zcat ${DATADIR}/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | \
    vcf-subset -c ${SAMPLE} | bgzip -c > ${PREFIX}.all.vcf.gz
tabix -f -p vcf ${PREFIX}.all.vcf.gz

# Get only 23andme
intersectBed -a ${PREFIX}.all.vcf.gz -b ${DATADIR}/23andme_snps_hg19.bed.gz -header | bgzip -c > \
    ${PREFIX}.subset.vcf.gz
tabix -f -p vcf ${PREFIX}.subset.vcf.gz
intersectBed -a ${PREFIX}.all.vcf.gz -b ${DATADIR}/23andme_snps_hg19.bed.gz -header -v | bgzip -c > \
    ${PREFIX}.heldout.vcf.gz
tabix -f -p vcf ${PREFIX}.heldout.vcf.gz

# Convert to gens file for IMPUTE2
vcf2impute_gen -vcf ${PREFIX}.subset.vcf.gz \
    -gen ${PREFIX}.subset.gen.gz 

# Convert heldout to gens file for IMPUTE2 to compare
vcf2impute_gen -vcf ${PREFIX}.heldout.vcf.gz \
    -gen ${PREFIX}.heldout.gen.gz 
