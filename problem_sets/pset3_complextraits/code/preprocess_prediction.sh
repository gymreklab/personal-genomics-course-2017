#!/bin/bash

KGPREFIX=/storage/resources/datasets/1000Genomes/phase3/
PREFIX=../data/ps3_pred

# Eye color
firstline=0
while IFS='' read -r line || [[ -n "$line" ]];
do
    chrom=$(echo $line | cut -f 1 -d' ')
    start=$(echo $line | cut -f 2 -d' ')
    end=$(echo $line | cut -f 3 -d' ')
    rsid=$(echo $line | cut -f 4 -d' ')
    if [ "$firstline" -eq "0" ]
    then
	args="-h"
    else
	args=""
    fi
    tabix $args ${KGPREFIX}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${chrom}:${start}-${start} | \
	awk -v"rsid=$rsid" '($3==rsid||$0~/#/)'
    firstline=1
done < ../data/eyecolor_snps.bed  | \
    vcf-subset -f -c ../data/CEUYRI.list | bgzip -c > ${PREFIX}_eyecolor.vcf.gz

# Eye color
firstline=0
while IFS='' read -r line || [[ -n "$line" ]];
do
    chrom=$(echo $line | cut -f 1 -d' ')
    start=$(echo $line | cut -f 2 -d' ')
    end=$(echo $line | cut -f 3 -d' ')
    rsid=$(echo $line | cut -f 4 -d' ')
    if [ "$firstline" -eq "0" ]
    then
	args="-h"
    else
	args=""
    fi
    tabix $args ${KGPREFIX}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${chrom}:${start}-${start} | \
	awk -v"rsid=$rsid" '($3==rsid||$0~/#/)'
    firstline=1
done < ../data/height_snps.bed  | \
    vcf-subset -f -c ../data/CEUYRI.list | bgzip -c > ${PREFIX}_height.vcf.gz

