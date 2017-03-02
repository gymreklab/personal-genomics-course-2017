#!/bin/bash

SUMMFILE=/oasis/projects/nsf/csd524/mgymrek/data/ps5/1kg_phase3_exome_ceu_all_vep.vcf.tab

#######################################
# Question 1: for sample NA06984
cat ${SUMMFILE} | grep "stop_gained\|frameshift" | awk '($7>0)' | cut -f 4 | sort | uniq
#######################################

#######################################
# Question 2A: rank each gene by number of samples with a nonsense/frameshift. Report top 10 genes. Where is MLL2?
opcols=$(seq 2 100 | awk '{print "sum " $0}')
cat ${SUMMFILE} | grep "stop_gained\|frameshift" | cut -f 1,2,3,5,6 --complement | \
    datamash -g 1 ${opcols} | ./count_greater_than_zero.py - 1 99 | \
    sort -k 2 -n -r

# Question 2B: rank each gene by number of samples with a missense
#######################################

#######################################
# Question 3
cat ${SUMMFILE} | awk '($6==0)' | grep "stop_gained\|frameshift" 

#######################################

#######################################
# Question 4: How many LoF MLL2 are in ExAC? What are their allele frequencies, where do they fall?
# see exac.broadinstitute.org
#######################################


#######################################
# Question 5: average MAF by category
for ann in  stop_gained synonymous_variant missense_variant
do
    echo $ann $(cat ${FINALOUTFILE}.tab | grep -v exacmaf | awk '($6>0)' \
	| awk -v"ann=$ann" '($5==ann)' | \
	cut -f 6 | grep -v "&" | sed 's/^.*://' | \
	sed 's/^$/0/' | datamash mean 1 )
done
