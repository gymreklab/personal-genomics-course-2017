#!/bin/bash

DATADIR=/oasis/projects/nsf/csd524/mgymrek/data/ps3/
PREFIX=${DATADIR}/ps3_gwas
USERPREFIX=/oasis/projects/nsf/csd524/$USER/ps3/results/ps3_gwas

# Run association test
plink \
    --bfile ${PREFIX} \
    --pheno ${PREFIX}.phen \
    --out ${USERPREFIX} \
    --linear \
    --allow-no-sex

# Change to more useable output format
cat ${USERPREFIX}.assoc.linear | sed -r 's/^\s+//g' | sed -r 's/\s+/\t/g' > ${USERPREFIX}.assoc.linear.tab
