#!/usr/bin/env python

"""
CSE291 - Personal Genomics for Bioinformaticians
Problem Set 2 - Ancestry: Principal components analysis

./ps3_manhattan.py \
  --assoc ps3_gwas.qassoc.tab \
  --out ps3_gwas.manhattan.png
"""

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--pcs", help="Tab-delim plink eigenvec results", type=str, required=True)
    parser.add_argument("--out", help="Output pdf or png file", type=str, required=True)
    args = parser.parse_args()
    
    # Load data
    eig = pd.read_csv(args.pcs, names=["sample", "family", "pc1", "pc2"], usecols=range(4), sep=" ")

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(eig.pc1, eig.pc2, color="blue", alpha=0.5)
    ax.set_xlabel("PC1", size=15)
    ax.set_ylabel("PC2", size=15)
    ax.set_xticklabels(ax.get_xticks(), size=12)
    ax.set_yticklabels(ax.get_yticks(), size=12)
    fig.savefig(args.out)

main()
