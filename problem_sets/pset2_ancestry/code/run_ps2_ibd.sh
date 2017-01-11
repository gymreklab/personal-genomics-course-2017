#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 1000
#SBATCH --get-user-env

DATADIR=/oasis/projects/nsf/csd524/mgymrek/data/ps2/
PREFIX=${DATADIR}/ps2_ibd
OUTDIR=/oasis/projects/nsf/csd524/$USER/ps2/results/

plink --bfile ${PREFIX}.lwk --genome --out ${OUTDIR}/lwk.ibd

