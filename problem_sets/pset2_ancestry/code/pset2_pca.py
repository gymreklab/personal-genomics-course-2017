#!/usr/bin/env python

"""
CSE291 - Personal Genomics for Bioinformaticians
Problem Set 2 - Ancestry: Principal components analysis

./pset2_pca.py \
  --ref-vcf ../data/ps2_reference_genotypes.vcf.gz \
  --ref-labels ../data/ps2_reference_labels.csv \
  --sample-vcf ../data/ps2_sample_genotypes.vcf.gz \
  --num-samples 100 \
  --num-snps 1000
"""

import argparse
import numpy as np
import sys
import vcf

def convert_vcf_to_matrix(vcffile, labelfile=None, \
                              numsamples=-1, numsnps=-1):
    """Function to convert VCF to matrix ready for PCA
    
    Args:
        vcffile (str): Path to input VCF file.
        labelfile (str): Path to population label CSV file. If
            not provided, use lable "Unknown".
        numsamples (int): Number of samples to use. -1 for all (default).
        numsnps (int): Number of SNPs to use. -1 for all (default).

    Returns:
        gtmatrix (np.array): Integer matrix of size numsamples x numsnps.
        samplelabels (list[str]): List of sample population labels.
    """
    return np.array([]), [] # TODO

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--ref-vcf", help="Reference panel VCF", type=str, required=True)
    parser.add_argument("--ref-labels", help="CSV file with sample, pop label for ref panel", type=str, required=True)
    parser.add_argument("--sample-vcf", help="VCF with samples of unknown ancestry", type=str, required=True)
    parser.add_argument("--num-samples", help="Number of samples to use from ref panel. -1 for all.", type=int, default=-1)
    parser.add_argument("--num-snps", help="Number of SNPs to include in analysis. -1 for all.", type=int, default=-1)
    parser.add_argument("--out", help="Prefix for output files.", type=str, default="ps2_pca")
    args = parser.parse_args()

    sys.stderr.write("[pset2_pca.py] Convert reference panel VCF to matrix...\n")
    ref_matrix, ref_labels = convert_vcf_to_matrix(args.ref_vcf, args.ref_labels, \
                                                       numsamples=args.num_samples, numsnps=args.num_snps)

    sys.stderr.write("[pset2_pca.py] Normalize SNP genotypes...\n") # TODO

    sys.stderr.write("[pset2_pca.py] Perform PCA...\n") # TODO

    sys.stderr.write("[pset2_pca.py] Plot reference panel...\n") # TODO

    sys.stderr.write("[pset2_pca.py] Analyze held out samples...\n") # TODO

    sys.stderr.write("[pset2_pca.py] Plot held out samples...\n") # TODO

    sys.stderr.write("[pset2_pca.py] Save plot and exit...\n") # TODO

if __name__ == "__main__":
    main()
