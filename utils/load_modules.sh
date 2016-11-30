#!/bin/bash

#################################################
# CSE291: Bioinformatics for personal genomics

# This script loads modules used in the course
#################################################

module load bcftools/1.3
module load bedtools/2.25.0
module load python/2.7.10
module load scipy/2.7
module load vcftools/0.1.14

export PATH=$PATH:/oasis/projects/nsf/csd524/mgymrek/bin
export PATH=$PATH:/home/mgymrek/workspace/personal-genomics-course-2017/problem_sets/pset2_ancestry/code/

