#!/usr/bin/env python

"""
CSE291 - Personal Genomics for Bioinformaticians
Problem Set 2 - Ancestry: Principal components analysis

./pset2_pca.py \
  --data-prefix ../data/ps2_pca \
  --labels ../data/ps2_reference_labels.csv \
  --num-samples 100 \
  --num-snps 1000
"""

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd # TODO maybe a way to load matrix without pandas
from sklearn.decomposition import PCA
import sys
import vcf

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

def normalize_snp_matrix(ref_matrix):
    """Function to normalize each SNP according to allele frequency
    
    Args:
        ref_matrix (np.array): Float matrix of size numsnps x numsamples

    Returns:
        ref_matrix_norm (np.array): Float matrix of normalized SNPs, size numsnps x numsamples
    """
    snp_p = np.mean(ref_matrix, axis=1, dtype=np.float64).reshape(ref_matrix.shape[0], 1)
    snp_se = np.sqrt(snp_p*(1-snp_p)).reshape(ref_matrix.shape[0], 1)
    ref_matrix_norm = (ref_matrix-snp_p)/snp_se
    return ref_matrix_norm

def perform_pca(ref_matrix_norm):
    """Perform PCA on normalized matrix

    Args:
        ref_matrix_norm (np.array): Float matrix of normalized SNPs, size numsnps x numsamples

    Returns:
        ref_matrix_projection (np.array): Float matrix of data projected onto PCs. size numsamples x numdim (2)
    """
    ndim = 2
    # Get GRM (covariance matrix)                                                                                                                                                             
    grm = ref_matrix_norm.transpose().dot(ref_matrix_norm)/ref_matrix_norm.shape[0]
    # Get Eigendecomposition - column i is ith eigenvector                                                                                                                                    
    evals, evecs = np.linalg.eig(grm)
    # sort eigenvalue in decreasing order
    idx = np.argsort(evals)[::-1]
    evecs = evecs[:,idx]
    # sort eigenvectors according to same index
    evals = evals[idx]
    # select the first n eigenvectors (n is desired dimension
    # of rescaled data array, or dims_rescaled_data)
    evecs = evecs[:, :ndim]
    return evecs

# Option 2: Using sklearn
#    pca = PCA(n_components=2)
#    pca.fit(ref_matrix_norm.transpose())
#    return pca.transform(ref_matrix_norm.transpose())
# Alternatively:
#    ndim = 2
#    # Get GRM (covariance matrix)                                                                                                                                                             
#    grm = ref_matrix_norm.dot(ref_matrix_norm.transpose())/ref_matrix_norm.shape[1]
#    # Get Eigendecomposition - column i is ith eigenvector                                                                                                                                    
#    evals, evecs = np.linalg.eig(grm)
#    # sort eigenvalue in decreasing order
#    idx = np.argsort(evals)[::-1]
#    evecs = evecs[:,idx]
#    # sort eigenvectors according to same index
#    evals = evals[idx]
#    # select the first n eigenvectors (n is desired dimension
#    # of rescaled data array, or dims_rescaled_data)
#    evecs = evecs[:, :ndim]
#    ref_matrix_projection = np.dot(evecs.transpose(), ref_matrix_norm).transpose()
#    return ref_matrix_projection


def plot_ref_panel(pc1, pc2, colors, samples, outprefix):
    """Function to plot top PCs of reference panel
    
    Args:
        pc1 (list[float]): First PC to plot
        pc2 (list[float]): Second PC to plot
        colors (list[str]): List of colors for each sample
        samples (list[str]): List of sample labels for each sample
        outprefix: prefix to name output file
    """
    pc1 = map(float, pc1)
    pc2 = map(float, pc2)
    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(pc1, pc2, color=colors, s=5, alpha=0.5)
    # Plot gray ones on top, print out their values
    idx = [i for i in range(len(colors)) if colors[i] == "gray"]
    for ind in idx:
        ax.scatter([pc1[ind]], [pc2[ind]], color="gray")
        ax.text(pc1[ind], pc2[ind], samples[ind])
        sys.stderr.write("\t".join(map(str, [samples[ind], pc1[ind], pc2[ind]]))+"\n")
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

def read_pca_data(dataprefix, labelfile, numsnps=-1):
    """Function to read matrix that was pre-calcluated
    Args:
        dataprefix (str): Prefix of files to load
        labelfile (str): Path to population label CSV file. If
            not provided, use lable "Unknown".
        numsnps (int): Number of SNPs to use. -1 for all (default).

    Returns:
        ref_matrix (np.array): Float matrix of size numsnps x numsamples
        ref_labels (list[str]): List of sample population labels.
        sample_labels (list[str]): List of sample IDs
    """
    ref_matrix = np.array(pd.read_csv(dataprefix+".genotypes.tab", sep="\t", header=False))
    sample_labels = map(lambda x: x.strip(), open(dataprefix + ".samples.txt", "r").readlines())
    # Get sample->pop dictionary
    if labelfile is not None:
        sample_to_pop = dict(map(lambda x: x.strip().split(","), open(labelfile, "r").readlines()))
    else: sample_to_pop = {}
    ref_labels = [sample_to_pop.get(s, "None") for s in sample_labels]
    if numsnps == -1: numsnps = ref_matrix.shape[0]
    return ref_matrix[0:numsnps,:], ref_labels, sample_labels
    
def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--num-samples", help="Number of samples to use from ref panel. -1 for all.", type=int, default=-1)
    parser.add_argument("--num-snps", help="Number of SNPs to include in analysis. -1 for all.", type=int, default=-1)
    parser.add_argument("--out", help="Prefix for output files.", type=str, default="ps2_pca")
    parser.add_argument("--data-prefix", help="Data prefix for precomputed genotype matrix", type=str, required=True)
    parser.add_argument("--labels", help="CSV file with sample, pop label for ref panel", type=str)
    args = parser.parse_args()

    sys.stderr.write("[pset2_pca.py] Read data for PCA...\n")
    ref_matrix, ref_labels, sample_labels = read_pca_data(args.data_prefix, args.labels, numsnps=args.num_snps)
    
    sys.stderr.write("[pset2_pca.py] Normalize SNP genotypes...\n")
    ref_matrix_norm = normalize_snp_matrix(ref_matrix)

    sys.stderr.write("[pset2_pca.py] Perform PCA...\n")
    ref_matrix_projection = perform_pca(ref_matrix_norm)

    sys.stderr.write("[pset2_pca.py] Plot reference panel...\n")
    plot_ref_panel(ref_matrix_projection[:,0], ref_matrix_projection[:, 1], \
                   [pop_to_color[pop] for pop in ref_labels], sample_labels, args.out)

if __name__ == "__main__":
    main()
