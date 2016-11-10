#!/usr/bin/env python

"""
CSE291 - Personal Genomics for Bioinformaticians 
Problem Set 2 - Ancestry: Principal components analysis                                                                                                                           
Convert vcf to tab to genotype matrix row by row

cat vcf-to-tab | ./get_vcf_matrix.py <OUTPREFIX>

"""

import numpy as np
import sys

MINMAF=0.01

def ProcessLine(line):
    """Convert row of ref, gt1, gt2, ... into genotype matrix row

    Args:
        line (str): input line

    Returns:
        newline (str): processed input line
    """
    items = line.strip().split()
    ref = items[2]
    genotypes = map(lambda x: ((1-x.split("/").count(ref)*0.5)), items[3:])
    if np.mean(genotypes) <= MINMAF or np.mean(genotypes) >= 1-MINMAF: return ""
    return "\t".join(map(str, genotypes))+"\n"

try:
    outprefix = sys.argv[1]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

# Write samples to file
samples = sys.stdin.readline().strip().split()[3:]
f = open(outprefix + ".samples.txt", "w")
for s in samples: f.write(s + "\n")
f.close()

# Write genotypes
f = open(outprefix + ".genotypes.tab", "w")
line = sys.stdin.readline()
while line != "":
    f.write(ProcessLine(line))
    f.flush()
    line = sys.stdin.readline()
f.close()

