#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 1000
#SBATCH --get-user-env

/oasis/projects/nsf/csd524/$USER/ps2/code/pset2_pca.py \
    --data-prefix /oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_pca \
    --labels /oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_reference_labels.csv \
    --out /oasis/projects/nsf/csd524/$USER/ps2/results/ps2_pca

#/oasis/projects/nsf/csd524/$USER/ps2/code/pset2_pca.py \
#    --data-prefix /oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_pca_23andme \
#    --labels /oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_reference_labels.csv \
#    --out /oasis/projects/nsf/csd524/$USER/ps2/results/ps2_pca_23andme
