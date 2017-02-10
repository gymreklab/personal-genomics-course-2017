#!/bin/bash

RESULTS=/oasis/projects/nsf/csd524/$USER/ps4/results
FAREF=/oasis/projects/nsf/csd524/mgymrek/data/ps4/GRCh38_full_analysis_set_plus_decoy_hla.fa
CRAMFILE=/oasis/projects/nsf/csd524/mgymrek/data/ps4/NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram

# Get pileup - for problem 2
samtools mpileup -r chr22 -f ${FAREF} ${CRAMFILE} > ${RESULTS}/NA12878.pileup

# Coverage histogram
samtools mpileup -f ${FAREF} -r chr1 ${CRAMFILE} | cut -f 4 | sort | uniq -c | tr -s ' ' '\t' | sed -e 's/^[ \t]*//' > ${RESULTS}/auto.cov
samtools mpileup -f ${FAREF} -r chrX ${CRAMFILE} | cut -f 4 | sort | uniq -c | tr -s ' ' '\t' | sed -e 's/^[ \t]*//' > ${RESULTS}/chrX.cov
samtools mpileup -f ${FAREF} -r chrY ${CRAMFILE} | cut -f 4 | sort | uniq -c | tr -s ' ' '\t' | sed -e 's/^[ \t]*//' > ${RESULTS}/chrY.cov

# Mate pair distances
samtools view ${CRAMFILE} chr1 | cut -f 9 | sort | uniq -c | tr -s ' ' '\t' | sed -e 's/^[ \t]*//' > ${RESULTS}/template_lens.txt

# Look at example reads
# samtools view $CRAMFILE chr1:949140-949150
# samtools tview ${CRAMFILE} ${FAREF}