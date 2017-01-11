#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 10
#SBATCH --get-user-env

for numsamples in 100 500 1000 2000 2500
do
    echo $numsamples
    time /oasis/projects/nsf/csd524/$USER/ps2/code/pset2_pca.py \
	--data-prefix /oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_pca \
	--labels /oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_reference_labels.csv \
	--out /oasis/projects/nsf/csd524/$USER/ps2/results/tmp \
	--num-samples $numsamples
done