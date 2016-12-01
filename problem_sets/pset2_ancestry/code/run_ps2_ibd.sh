#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 10
#SBATCH --get-user-env

#java -jar -Xmx512m /home/mgymrek/sources/beagle.r1399.jar \
#    gt=/oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_ibd.vcf.gz \
#    out=/oasis/projects/nsf/csd524/$USER/ps2/results/ps2_ibd \
#    ibd=true impute=false

#cat /oasis/projects/nsf/csd524/$USER/ps2/results/ps2_ibd.ibd | \
#    java -jar -Xmx512m /home/mgymrek/sources/ibdmerge.jar > \
#    /oasis/projects/nsf/csd524/$USER/ps2/results/ps2_ibd.ibd.merged

plink --tfile /oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_ibd \
    --genome \
    --out /oasis/projects/nsf/csd524/$USER/ps2/results/ps2_ibd_plink