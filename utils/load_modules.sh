#!/bin/bash

#################################################
# CSE291: Bioinformatics for personal genomics

# This script loads modules used in the course
#################################################

module load bcftools/1.3
module load bedtools/2.25.0
module load plink/1.9
module load python/2.7.10
module load scipy/2.7
module load vcftools/0.1.14
module load bwa/0.7.13

export PATH=$PATH:/oasis/projects/nsf/csd524/mgymrek/bin
export PATH=$PATH:/home/mgymrek/workspace/personal-genomics-course-2017/problem_sets/pset2_ancestry/code/

# For cram to bam
export REF_PATH=/oasis/projects/nsf/csd524/mgymrek/data/ps4/cramcache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
export REF_CACHE=/oasis/projects/nsf/csd524/mgymrek/data/ps4/cramcache/%2s/%2s/%s