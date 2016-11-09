#!/usr/bin/env python

"""
CSE291 - Personal Genomics for Bioinformaticians
Problem Set 2 - Ancestry: Principal components analysis

./pset2_pca.py \
  --vcf ../data/ps2_reference_genotypes.vcf.gz \
  --labels ../data/ps2_reference_labels.csv \
  --num-samples 100 \
  --num-snps 1000
"""

import argparse
import cPickle
import matplotlib.pyplot as plt
import numpy as np
import sys
import vcf

MINMAF = 0.1

pop_to_color = {
    "ACB": "blue",
    "ASW": "blue",
    "BEB": "blue",
    "CDX": "green",
    "CEU": "yellow",
    "CHB": "green",
    "CHS": "green",
    "CLM": "purple",
    "ESN": "yellow",
    "FIN": "yellow",
    "GBR": "yellow",
    "GIH": "orange",
    "GWD": "orange",
    "IBS": "yellow",
    "ITU": "yellow",
    "JPT": "green",
    "KHV": "red",
    "LWK": "red",
    "MSL": "red",
    "MXL": "purple",
    "PEL": "purple",
    "PJL": "purple",
    "PUR": "purple",
    "STU": "purple",
    "TSI": "yellow",
    "YRI": "red",
    "None": "gray"
}

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
        gtmatrix (np.array): Float matrix of size numsamples x numsnps.
                 0, 0.5, 1 for AA, AB, BB
        samplelabels (list[str]): List of sample population labels.
    """
    # Get sample->pop dictionary
    if labelfile is not None:
        sample_to_pop = dict(map(lambda x: x.strip().split(","), open(labelfile, "r").readlines()))
    else: sample_to_pop = {}
    # Read SNP genotype matrix - TODO preallocate space for this to save time?
    reader = vcf.Reader(open(vcffile, "rb"))
    if numsamples == -1 or numsamples > len(reader.samples): numsamples = len(reader.samples)
    snp_gts = []
    for record in reader:
        if len(record.alleles) > 2: continue
        if record.aaf[0] <= MINMAF or record.aaf[0] >= (1-MINMAF): continue
        snpdata = [sum((map(float, sample.gt_alleles)))*0.5 for sample in record]
        snp_gts.append(snpdata)
        if len(snp_gts) >= numsnps and numsnps != -1: break
        if len(snp_gts)%100 == 1: sys.stderr.write("[pset2_pca.py]    Processed %s SNPs...\n"%len(snp_gts))
    # Return results
    return np.array(snp_gts)[:,0:numsamples], \
        map(lambda x: sample_to_pop.get(x, "None"), reader.samples[0:numsamples])

def normalize_snp_matrix(ref_matrix):
    """Function to normalize each SNP according to allele frequency
    
    Args:
        ref_matrix (np.array): Float matrix of size numsamples x numsnps

    Returns:
        ref_matrix_norm (np.array): Float matrix of normalized SNPs, size numsamples x numsnps
        snp_p (np.array): Float array of allele frequencies. size numsnps
    """
    snp_p = np.mean(ref_matrix, axis=1, dtype=np.float64).reshape(ref_matrix.shape[0], 1)
    snp_se = np.sqrt(snp_p*(1-snp_p)).reshape(ref_matrix.shape[0], 1)
    ref_matrix_norm = (ref_matrix-snp_p)/snp_se
    return ref_matrix_norm, snp_p

def perform_pca(ref_matrix_norm):
    """Perform PCA on normalized matrix

    Args:
        ref_matrix_norm (np.array): Float matrix of normalized SNPs, size numsamples x numsnps

    Returns:
        ref_matrix_projection (np.array): Float matrix of data projected onto PCs. size numsamples x numsnps
        eigvecs (np.array): Matrix of eigenvectors. size numsamples x numsamples
    """
    # TODO are these sorted?
    # Get GRM (covariance matrix)
    grm = ref_matrix_norm.transpose().dot(ref_matrix_norm)/ref_matrix_norm.shape[0]
    # Get Eigendecomposition
    eigvals, eigvecs = np.linalg.eig(grm)
    # Project data onto PCs
    ref_matrix_projection = eigvecs.transpose().dot(ref_matrix_norm.transpose())
    return ref_matrix_projection, eigvecs

def plot_ref_panel(pc1, pc2, colors, outprefix):
    """Function to plot top PCs of reference panel
    
    Args:
        pc1 (list[float]): First PC to plot
        pc2 (list[float]): Second PC to plot
        colors (list[str]): List of colors for each sample
        outprefix: prefix to name output file
    """
    pc1 = map(float, pc1)
    pc2 = map(float, pc2)
    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(pc1, pc2, color=colors)
    # Now make it look a little prettier since I can't stand matplotlib defaults...
    ax.set_xlabel("PC 1", size=15)
    ax.set_ylabel("PC 2", size=15)
    ax.set_xticklabels(ax.get_xticks(), size=12)
    ax.set_yticklabels(ax.get_yticks(), size=12)
    ax.spines["top"].set_visible(False);
    ax.spines["right"].set_visible(False);
    ax.get_xaxis().tick_bottom();
    ax.get_yaxis().tick_left();
    fig.savefig(outprefix + ".pdf")

def read_pca_data(outprefix, numsamples=-1):
    """Function to read matrix that was pre-calcluated
    Args:
        outprefix (str): Prefix of files to load
        numsamples (int): Number of samples to use. -1 for all (default).

    Returns:
        ref_matrix (np.array): Float matrix of size numsamples x numsnps
        ref_labels (list[str]): List of sample population labels.
    """
    ref_matrix = cPickle.load(open(outprefix+".ref_matrix.pkl", "r"))
    ref_labels = cPickle.load(open(outprefix+".ref_labels.pkl", "r"))
    if numsamples != -1:
        ref_matrix = ref_matrix[:,0:numsamples]
        ref_labels = ref_labels[0:numsamples]
    return ref_matrix, ref_labels
    
def write_pca_data(ref_matrix, ref_labels, outprefix):
    """Function to write matrix so we don't have to recalculate each time.

    Args:
        ref_matrix (np.array): Float matrix of size numsamples x numsnps
        ref_labels (list[str]): List of sample population labels. 
        outprefix (str): Output prefix for files
    """
    cPickle.dump(ref_matrix, open(outprefix+".ref_matrix.pkl", "w"))
    cPickle.dump(ref_labels, open(outprefix+".ref_labels.pkl", "w"))
    
def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="Reference panel VCF", type=str, required=True)
    parser.add_argument("--labels", help="CSV file with sample, pop label for ref panel", type=str, required=True)
    parser.add_argument("--num-samples", help="Number of samples to use from ref panel. -1 for all.", type=int, default=-1)
    parser.add_argument("--num-snps", help="Number of SNPs to include in analysis. -1 for all.", type=int, default=-1)
    parser.add_argument("--out", help="Prefix for output files.", type=str, default="ps2_pca")
    parser.add_argument("--load-from-pkl", help="Load matrix from pickle (for debugging)", action="store_true")
    args = parser.parse_args()

    sys.stderr.write("[pset2_pca.py] Convert reference panel VCF to matrix...\n")
    if args.load_from_pkl:
        ref_matrix, ref_labels = read_pca_data(args.out, numsamples=args.num_samples)
    else:
        ref_matrix, ref_labels = convert_vcf_to_matrix(args.vcf, args.labels, \
                                                       numsamples=args.num_samples, numsnps=args.num_snps)
        write_pca_data(ref_matrix, ref_labels, args.out)
    
    sys.stderr.write("[pset2_pca.py] Normalize SNP genotypes...\n")
    ref_matrix_norm, allele_freqs = normalize_snp_matrix(ref_matrix)
    
    sys.stderr.write("[pset2_pca.py] Perform PCA...\n")
    ref_matrix_projection, eigvecs = perform_pca(ref_matrix_norm)

    sys.stderr.write("[pset2_pca.py] Plot reference panel...\n")
    plot_ref_panel(ref_matrix_projection[:,0], ref_matrix_projection[:, 1], \
                   [pop_to_color[pop] for pop in ref_labels], args.out)

if __name__ == "__main__":
    main()
