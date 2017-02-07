#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 1000
#SBATCH --get-user-env

DATADIR=/oasis/projects/nsf/csd524/mgymrek/data/ps4/
INDEX=${DATADIR}/hs37d5.fa
FQ1=${DATADIR}/SRR622457_1.fastq.gz
FQ2=${DATADIR}/SRR622457_2.fastq.gz

#bwa mem -R '@RG\tID:NA12878\tSM:NA12878' ${INDEX} ${FQ1} ${FQ2} > ${DATADIR}/NA12878.sam

# Convert to BAM, note sam is truncated but we'll just use a subset?
cat ${DATADIR}/NA12878.sam | grep -v "SRR622457\.68907082" | \
    samtools view -bS - > ${DATADIR}/NA12878.bam
