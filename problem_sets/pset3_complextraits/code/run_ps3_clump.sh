#!/bin/bash

DATADIR=/oasis/projects/nsf/csd524/mgymrek/data/ps3/
PREFIX=${DATADIR}/ps3_gwas
USERPREFIX=/oasis/projects/nsf/csd524/$USER/ps3/results/ps3_gwas

# Get independent signals
plink \
    --bfile ${PREFIX} \
    --clump ${USERPREFIX}.assoc.linear --clump-field P \
    --clump-p1 0.0001 \
    --clump-p2 0.01 \
    --clump-r2 0.5 \
    --clump-kb 250 \
    --out ${USERPREFIX}