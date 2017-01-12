#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 1000
#SBATCH --get-user-env

chrom=2
numsnps=25
DATADIR=/oasis/projects/nsf/csd524/mgymrek/data/ps3/
KGVCF=${DATADIR}/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

PREFIX=${DATADIR}/ps3_gwas

# Subset to get CEU, TSI
zcat ${KGVCF} | \
    grep -v "VT=SV" | \
    intersectBed -header -a stdin -b ${DATADIR}/../ps2/23andme_snps_hg19.bed.gz > \
    ${PREFIX}.subset.vcf
cat ${PREFIX}.subset.vcf | \
    vcf-subset -f -c ../data/CEUTSI.list | \
    uniq | bgzip -c > ${PREFIX}.vcf.gz
tabix -f -p vcf ${PREFIX}.vcf.gz

# Convert to plink
plink --vcf ${PREFIX}.vcf.gz \
    --make-bed \
    --out ${PREFIX}
plink --bfile ${PREFIX} --recode12 --out ${PREFIX}

# Get list of causal SNPs - randomly select from VCF
zcat ${PREFIX}.vcf.gz | grep -v "^#" | cut -f 3 | shuf | head -n ${numsnps} > ../data/causal.snplist

# Simulate phenotype using GCTA
gcta64 --bfile ${PREFIX} --simu-qt --simu-causal-loci ../data/causal.snplist --simu-hsq 0.8 --simu-rep 1 --out ${PREFIX}.raw

# Add covariate for population
./add_phenotype_covar.py ../data/CEUTSI.labels ${PREFIX}.raw.phen 0.5 > ${PREFIX}.phen
