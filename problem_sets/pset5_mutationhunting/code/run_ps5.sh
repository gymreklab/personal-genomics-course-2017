#!/bin/bash

SUMMFILE=/oasis/projects/nsf/csd524/mgymrek/data/ps5/1kg_phase3_exome_ceu_all_vep.vcf.tab

# Question 1: for sample NA06984 Answer: 95
cat ${SUMMFILE} | grep "stop_gained\|frameshift" | awk '($7>0)' | cut -f 4 | sort | uniq | wc -l

# Question 2A: rank each gene by number of samples with a nonsense/frameshift
opcols=$(seq 2 100 | awk '{print "sum " $0}')
cat ${SUMMFILE} | grep "stop_gained\|frameshift" | cut -f 1,2,3,5,6 --complement | \
    datamash -g 1 ${opcols} | ./count_greater_than_zero.py - 1 99 | \
    sort -k 2 -n -r
# Report top 10 genes. Where is MLL2?
# head -n 10
# ZNF598
# ZNF516
# ZNF274
# ZFPM1
# UNC93B1
# SYNM
# SCARF2
# REC8
# PRDM15
# PRB3
# grep -n KMT2D (# 64)

# Question 2B: rank each gene by number of samples with a missense
cat ${SUMMFILE} | grep "missense" | cut -f 1,2,3,5,6 --complement | \
    datamash -g 1 ${opcols} | ./count_greater_than_zero.py - 1 99 | \
    sort -k 2 -n -r

# Report top 10 genes. Where is MLL2?
#ZSCAN20
#ZNF787
#ZNF778
#ZNF708
#ZNF681
#ZNF675
#ZNF668
#ZNF667
#ZNF665
#ZNF607
# (# 1445)

# Repeat above, filter out anything in ExaC (use exac_maf)
cat ${SUMMFILE} | awk '($6==0)' | grep "stop_gained\|frameshift" | cut -f 1,2,3,5,6 --complement | \
    datamash -g 1 ${opcols} | ./count_greater_than_zero.py - 1 99 | \
    sort -k 2 -n -r | head -n 10
# RP11-1407O15.2
# KMT2D (#2)
# GSTT2
# PRAMEF4
# LILRB3
# GNRH2
# FOLR3
# TRAV26-1
# SLC22A9
# KIR2DS4

cat ${SUMMFILE} | awk '($6==0)' | grep "missense" | cut -f 1,2,3,5,6 --complement | \
    datamash -g 1 ${opcols} | ./count_greater_than_zero.py - 1 99 | \
    sort -k 2 -n -r | head -n 10
# LILRB3
# NBPF10
# RHD
# NBPF14
# KIR3DL1
# IGHV3-30
# KRT34
# MKI67
# KMT2D
# FCGBP

# Question 4: How many LoF in MLL2 in ExAC? (use the browser). what are allele freqs? where do they fall?
# 16. all super rare (allele count 1 for all except 3 and 4)

# Question 5: what are average allele freqs for nonsense, missense, synonymous, intron
for ann in  stop_gained synonymous_variant missense_variant frameshift_variant inframe_deletion
do
    echo $ann $(cat ${FINALOUTFILE}.tab | grep -v exacmaf | awk -v"ann=$ann" '($5==ann)' | cut -f 6 | grep -v "&" | sed 's/^.*://' | sed 's/^$/0/' | datamash mean 1 )
done
