#!/bin/bash

DATADIR=/oasis/projects/nsf/csd524/mgymrek/data/ps3/
RESULTSDIR=/oasis/projects/nsf/csd524/mgymrek/ps3/results/

# Output per-individual eyecolor
./pset3_predict.py \
    --genotypes ${DATADIR}/ps3_pred_eyecolor.vcf.gz \
    --snps ${DATADIR}/eyecolor_snps_irisplex.bed > \
    ${RESULTSDIR}/eyecolor_predictions.tab