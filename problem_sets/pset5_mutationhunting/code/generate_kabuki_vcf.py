#!/usr/bin/env python

"""
Generate fake VCF data for kabuki mutations
Usage:
./generate_kabuki_vcf.py <mutfile> <numsamples>
"""

import sys

try:
    mutfile = sys.argv[1]
    numsamples = int(sys.argv[2])
except:
    sys.stderr.write(__doc__)

mutcount = 0
with open(mutfile, "r") as f:
    for line in f:
        chrom, start, end, ref, alt, annotation = line.strip().split()
        vcf_line = "\t".join([chrom.strip("chr"), start, ".", ref, alt, "100", "PASS", "VT=SNP", "GT"])
        gts = []
        for i in range(numsamples):
            if i%5 == mutcount:
                gts.append("0|1")
            else: gts.append("0|0")
        mutcount += 1
        sys.stdout.write(vcf_line + "\t" + "\t".join(gts)+"\n")
                             
