#!/usr/bin/env python

"""
Preprocess VEP VCF to a file with columns:
chrom, start, rsid, genesymbol, consequence, exacmaf, gtcols of 0,1,2
Usage:
./preprocess_vepmatrix.py <vepvcf> > <veptab>
"""

import vcf
import sys
import gzip

try:
    invcf = sys.argv[1]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

reader = vcf.Reader(open(invcf, "rb"))
header = ["chrom", "start", "rsid", "gene", "consequence", "exacmaf"] + list(reader.samples)
sys.stdout.write("\t".join(header)+"\n")

with gzip.open(invcf, "r") as f:
    for line in f:
        if line.startswith("#"): continue
        items = line.strip().split()
        chrom = items[0]
        start = items[1]
        rsid = items[2]
        info = items[7].split(";")
        csq = []
        for item in info:
            if item.startswith("CSQ="):
                csq = item.split("|")
        if len(csq) == 0: continue
        gene = csq[3]
        consequence = csq[1]
        exac_maf = csq[45]
        pitems = [chrom, start, rsid, gene, consequence, exac_maf]
        for gt in items[9:]:
            try:
                pitems.append(sum(map(int, gt.split("|"))))
            except: pitems.append(0)
        sys.stdout.write("\t".join(map(str, pitems))+"\n") 
"""
for record in reader:
    chrom = record.CHROM
    start = record.POS
    rsid = record.ID
    csq = record.INFO["CSQ"]
    gene = csq[0].split("|")[3]
    consequence = csq[0].split("|")[1]
    exac_maf = csq[0].split("|")[45]
    items = [chrom, start, rsid, gene, consequence, exac_maf]
    for sample in record:
        items.append(sum(map(int, sample.gt_alleles)))
    sys.stdout.write("\t".join(map(str, items))+"\n")
"""
