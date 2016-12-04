#!/bin/bash

#SBATCH -A csd524
#SBATCH -p compute
#SBATCH -t 1000
#SBATCH --get-user-env

REFDIR=/oasis/projects/nsf/csd524/mgymrek/data/ps2/imputation/1000GP_Phase3
TMPLOC=/oasis/scratch/comet/$USER/temp_project
tmpdir=$(mktemp -d -p ${TMPLOC})
chrom=16

# Make separate CEU, YRI, CEU+YRI, and all ref panels
cat ${REFDIR}/1000GP_Phase3.sample | grep -n CEU | cut -f 1 -d':' > ${tmpdir}/ceu.ind
cat ${REFDIR}/1000GP_Phase3.sample | grep -n YRI | cut -f 1 -d':' > ${tmpdir}/yri.ind
cat ${REFDIR}/1000GP_Phase3.sample | grep -n "CEU\|YRI" | cut -f 1 -d':' > ${tmpdir}/ceu_yri.ind
cat ${REFDIR}/1000GP_Phase3.sample | grep -v -n ASW | cut -f 1 -d':' > ${tmpdir}/not_asw.ind

# Subset hap files
for ind in ceu yri ceu_yri not_asw
do
    indfile=${tmpdir}/${ind}.ind
    echo $ind
    cols=$(cat ${indfile} | awk '{print 2*$1-1 "\n" 2*$1}' | datamash transpose -t ',')
    zcat ${REFDIR}/1000GP_Phase3_chr${chrom}.hap.gz | \
	cut -f ${cols} -d' ' > ${REFDIR}/1000GP_Phase3_${ind}_chr${chrom}.hap
    gzip ${REFDIR}/1000GP_Phase3_${ind}_chr${chrom}.hap
done



