#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 1000
#SBATCH --get-user-env

RESULTS=/oasis/projects/nsf/csd524/$USER/ps4/results
FAREF=/oasis/projects/nsf/csd524/mgymrek/data/ps4/GRCh38_full_analysis_set_plus_decoy_hla.fa
CRAMFILE=/oasis/projects/nsf/csd524/mgymrek/data/ps4/NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram

# Run SNP caller
samtools mpileup -r chr22 -f ${FAREF} ${CRAMFILE} | grep -v "+" | grep -v "-" | \
    ./ps4_snpcaller.py - ${RESULTS}/NA12878_chr22_snps.tab
