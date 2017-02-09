#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 2000
#SBATCH --get-user-env

CRAMFILE=/oasis/projects/nsf/csd524/mgymrek/data/ps4/NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram
BAMFILE=/oasis/projects/nsf/csd524/mgymrek/data/ps4/NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.bam

samtools view -b ${CRAMFILE} > ${BAMFILE}
samtools index ${BAMFILE}