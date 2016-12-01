#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 1000
#SBATCH --get-user-env

INVCF=/oasis/projects/nsf/csd524/mgymrek/data/ps2/ALL.chr16.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz
SUBVCF=/oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_ibd.vcf.gz
PREFIX=/oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_ibd

#	vcf-subset -c NA19713,NA19985,NA19625,NA20414,NA19660,NA19685 \

# Get subset of samples
for chrom in $(seq 1 22)
do
    zcat /oasis/projects/nsf/csd524/mgymrek/data/ps2/ALL.chr${chrom}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz | grep -v "VT=SV" | \
	intersectBed -header -a stdin -b 23andme_snps_hg19.bed.gz \
	 > /oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_ibd_${chrom}.vcf
    bgzip -f /oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_ibd_${chrom}.vcf
    tabix -f -p vcf /oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_ibd_${chrom}.vcf.gz
done

# concat chrom VCFs
vcf-concat $(ls /oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_ibd_*.vcf.gz) | bgzip -c > ${SUBVCF}
tabix -f -p vcf ${SUBVCF}

# Convert to plink format
vcftools --gzvcf ${SUBVCF} --plink-tped --out ${PREFIX}
