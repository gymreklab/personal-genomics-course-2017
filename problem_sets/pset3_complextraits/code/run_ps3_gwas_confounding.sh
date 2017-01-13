#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 1000
#SBATCH --get-user-env

DATADIR=/oasis/projects/nsf/csd524/mgymrek/data/ps3/
PREFIX=${DATADIR}/ps3_gwas
USERPREFIX=/oasis/projects/nsf/csd524/$USER/ps3/results/ps3_gwas

# Use plink to get PCs
plink \
    --bfile ${PREFIX} \
    --pca 10 \
    --out ${USERPREFIX}

# Now do the GWAS with the covariate file
plink \
    --bfile ${PREFIX} \
    --pheno ${PREFIX}.phen \
    --out ${USERPREFIX}.withcovar \
    --covar ${USERPREFIX}.eigenvec \
    --linear \
    --allow-no-sex

# Change to more useable output format
cat ${USERPREFIX}.withcovar.assoc.linear | sed -r 's/^\s+//g' | sed -r 's/\s+/\t/g' | \
    grep "ADD\|CHR" > ${USERPREFIX}.withcovar.assoc.linear.tab
