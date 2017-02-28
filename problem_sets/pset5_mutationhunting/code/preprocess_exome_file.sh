#!/bin/bash

KGDIR=/oasis/projects/nsf/csd524/mgymrek/data/ps5/
DATADIR=/oasis/projects/nsf/csd524/mgymrek/data/ps5/

zcat ${KGDIR}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | head -n 10000 | \
    grep "^#" > ${DATADIR}/1kg_phase3_exome.vcf

for chrom in $(seq 1 22)
do
    intersectBed -a ${KGDIR}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
	-b ${KGDIR}/20130108.exome.targets_nochr.bed | grep -v "VT=SV"
done >> ${DATADIR}/1kg_phase3_exome.vcf

bgzip ${DATADIR}/1kg_phase3_exome.vcf
tabix -p vcf ${DATADIR}/1kg_phase3_exome.vcf.gz

# Subset individuals
zcat ${DATADIR}/1kg_phase3_exome.vcf.gz | \
    vcf-subset -f -c ../data/CEU.samples | bgzip -c > \
    ${DATADIR}/1kg_phase3_exome_ceu.vcf.gz
tabix -p vcf ${DATADIR}/1kg_phase3_exome_ceu.vcf.gz

# Add disease mutations (5 possible mutations. sample gets mutation corresponding to colnum%5
./generate_kabuki_vcf.py ../data/kabuki_mutations.bed 99 > kabuki_mutations.vcf
zcat ${DATADIR}/1kg_phase3_exome_ceu.vcf.gz | cat - kabuki_mutations.vcf >> ${DATADIR}/1kg_phase3_exome_ceu_v2.vcf

# Sort and index new file
cat ${DATADIR}/1kg_phase3_exome_ceu_v2.vcf | vcf-sort | bgzip -c > ${DATADIR}/1kg_phase3_exome_ceu_v2.vcf.gz
tabix -p vcf ${DATADIR}/1kg_phase3_exome_ceu_v2.vcf.gz