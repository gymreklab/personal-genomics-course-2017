#!/usr/bin/env python

"""
CSE291 - Personal Genomics for Bioinformaticians
Problem Set 2 - Ancestry: Principal components analysis

./ps3_manhattan.py \
  --assoc ps3_gwas.assoc.linear.tab \
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
    parser.add_argument("--assoc", help="Tab-delim plink association results", type=str, required=True)
    parser.add_argument("--out", help="Output pdf or png file", type=str, required=True)
    args = parser.parse_args()
    
    # Load data
    assoc = pd.read_csv(args.assoc, sep="\t")

    # Remove NAs
    assoc = assoc[~np.isnan(assoc["P"])]

    # Sort by position (should be already)
    assoc = assoc.sort("BP")

    # Plot
    fig = plt.figure()
    fig.set_size_inches((10, 5))
    ax = fig.add_subplot(111)
    pvals = -1*assoc["P"].apply(lambda x: np.log10(x))
    ax.scatter(assoc["BP"], pvals, color="blue", alpha=0.5)
    ax.set_xlabel("Position (chr2)", size=15)
    ax.set_ylabel("-log10 P", size=15)
    ax.set_xlim(left=0, right=max(assoc["BP"]))
    ax.set_ylim(bottom=0, top=max(pvals))
    ax.set_xticklabels(ax.get_xticks(), size=12)
    ax.set_yticklabels(ax.get_yticks(), size=12)
    fig.savefig(args.out)

main()
