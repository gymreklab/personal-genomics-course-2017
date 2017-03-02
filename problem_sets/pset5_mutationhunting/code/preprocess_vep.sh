#!/bin/bash

VCF=/oasis/projects/nsf/csd524/mgymrek/data/ps5/1kg_phase3_exome_ceu.vcf.gz
OUTFILE=/oasis/projects/nsf/csd524/mgymrek/data/ps5/1kg_phase3_exome_ceu_vep_1kg.vcf
VEPPATH=/home/mgymrek/sources/ensembl-tools-release-87/scripts/variant_effect_predictor/
KVCF=kabuki_mutations.vcf
KOUTFILE=/oasis/projects/nsf/csd524/mgymrek/data/ps5/1kg_phase3_exome_ceu_kabuki_vep.vcf
FINALOUTFILE=/oasis/projects/nsf/csd524/mgymrek/data/ps5/1kg_phase3_exome_ceu_all_vep.vcf

perl ${VEPPATH}/variant_effect_predictor.pl \
    --everything \
    --offline \
    --vcf --force_overwrite \
    -i ${VCF} \
    -o ${OUTFILE} \
    --chr 17,18,19,20,21,22

# For Kabuki mutations separately
perl ${VEPPATH}/variant_effect_predictor.pl \
    --everything --offline --vcf --force_overwrite \
    -i ${KVCF} \
    -o ${KOUTFILE}

head -n 10000 ${OUTFILE} | grep "^#" > ${FINALOUTFILE}
cat ${OUTFILE} ${KOUTFILE} | grep -v "^#" >> ${FINALOUTFILE}
cat ${FINALOUTFILE} | vcf-sort | bgzip -c > ${FINALOUTFILE}.gz
tabix -p vcf ${FINALOUTFILE}.gz

# Convert output to more easily parseable tab format
./preprocess_vepmatrix.py ${FINALOUTFILE}.gz | uniq > ${FINALOUTFILE}.tab

# Now groupby each mutation type, count number in each individual
cat ${FINALOUTFILE}.tab | head -n 1 | cut -f 1-3,6 --complement > ${FINALOUTFILE}.summ.tab
colops=$(seq 3 101 | awk '{print "sum " $0}')
cat ${FINALOUTFILE}.tab | grep -v exacmaf | cut -f 1-3,6 --complement | \
    sort -k1,1 -k2,2 | datamash -g 1,2 ${colops} >> ${FINALOUTFILE}.summ.tab
