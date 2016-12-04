#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 1000
#SBATCH --get-user-env

DATADIR=/oasis/projects/nsf/csd524/mgymrek/data/ps2/
OUTDIR=/oasis/projects/nsf/csd524/mgymrek/ps2/results/
chrom=16

# Run with 4 different reference panels
for refpanel in ceu yri ceu_yri not_asw
do
    impute2 \
	-m ${DATADIR}/imputation/1000GP_Phase3/genetic_map_chr${chrom}_combined_b37.txt \
	-h ${DATADIR}/imputation/1000GP_Phase3/1000GP_Phase3_${refpanel}_chr${chrom}.hap.gz \
	-l ${DATADIR}/imputation/1000GP_Phase3/1000GP_Phase3_chr${chrom}.legend.gz \
	-g ${DATADIR}/ps2_impute.subset.gen.gz \
	-int 5e6 10e6 \
	-Ne 20000 \
	-o ${OUTDIR}/ps2_impute.phased.${refpanel}.impute2
done