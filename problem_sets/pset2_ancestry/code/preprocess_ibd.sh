#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 1000
#SBATCH --get-user-env

DATADIR=/oasis/projects/nsf/csd524/mgymrek/data/ps2/
PREFIX=${DATADIR}/ps2_ibd

# Get VCF
for chrom in $(seq 1 22)
do
    zcat ${DATADIR}/ALL.chr${chrom}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz | \
	grep -v "VT=SV" | \
	intersectBed -header -a stdin -b 23andme_snps_hg19.bed.gz \
	 > ${PREFIX}_${chrom}.vcf
    bgzip -f ${PREFIX}_${chrom}.vcf
    tabix -f -p vcf ${PREFIX}_${chrom}.vcf.gz
done

# concat chrom VCFs
vcf-concat $(ls ${PREFIX}_*.vcf.gz) | bgzip -c > ${PRERIX}.vcf.gz
tabix -f -p vcf ${PREFIX}.vcf.gz

# Subset
zcat ${PREFIX}.vcf.gz | vcf-subset -f -c LWK.list | uniq | bgzip -c > ${PREFIX}.lwk.vcf.gz
tabix -f -p vcf ${PREFIX}.lwk.vcf.gz
plink --vcf ${PREFIX}.lwk.vcf.gz \
    --cm-map ${DATADIR}/imputation/1000GP_Phase3/genetic_map_chr@_combined_b37.txt \
    --make-bed \
    --out ${PREFIX}.lwk

# Remove chrom files
rm ${PREFIX}_*.vcf.gz
rm ${PREFIX}_*.vcf.gz.tbi
