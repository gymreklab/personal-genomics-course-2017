#!/usr/bin/env python

"""
CSE291 - Personal Genomics for Bioinformaticians
Problem Set 3 - Complex traits

./pset3_predict.py \
  --genotypes VCFFILE
  --snps SNPBEDFILE

computes pi1, pi2, pi3 using IrisPlex predictions
"""

alpha1 = 3.94
alpha2 = 0.65

import argparse
import math
import pandas as pd
import sys
import vcf

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--genotypes", help="VCF file with sample genotypes", required=True, type=str)
    parser.add_argument("--snps", help="BED file with SNPs for prediction", required=True, type=str)
    args = parser.parse_args()

    # Load SNPs
    snps = pd.read_csv(args.snps, sep=r"\s*", names=["chrom","start","end","rsid","effect_allele","b1","b2"])

    # Keep scores for each sample in the VCF file
    vcfreader = vcf.Reader(open(args.genotypes, "rb"))
    sample_to_score1 = dict(zip(vcfreader.samples, [0]*len(vcfreader.samples)))
    sample_to_score2 = dict(zip(vcfreader.samples, [0]*len(vcfreader.samples)))

    # Go through each SNP
    for i in range(snps.shape[0]):
        chrom = snps["chrom"].values[i]
        start = int(snps["start"].values[i])
        eff_size1 = float(snps["b1"].values[i])
        eff_size2 = float(snps["b2"].values[i])
        eff_allele = snps["effect_allele"].values[i]
        rsid = snps["rsid"].values[i]
        matches = vcfreader.fetch("%s:%s"%(chrom, start))
        for snp in matches: # should only be one
            alleles = map(str, snp.alleles)
            for sample in snp:
                num_eff_alleles = sum(map(lambda x: alleles[int(x)]==eff_allele, sample.gt_alleles))
                print eff_size1*num_eff_alleles, eff_size2*num_eff_alleles
                sample_to_score1[sample.sample] += eff_size1*num_eff_alleles
                sample_to_score2[sample.sample] += eff_size2*num_eff_alleles

    for sample in sample_to_score1:
        x1 = math.exp(alpha1 + sample_to_score1[sample])
        x2 = math.exp(alpha2 + sample_to_score2[sample])
        pi1 = x1/(1+x1+x2)
        pi2 = x2/(1+x1+x2)
        pi3 = 1-pi1-pi2
        sys.stdout.write("\t".join(map(str, [sample, pi1, pi2, pi3]))+"\n")

if __name__ == "__main__":
    main()
