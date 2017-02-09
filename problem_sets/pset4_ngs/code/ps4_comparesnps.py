#!/usr/bin/env python

"""
CSE291 SNP calling execise - problem set 4

Usage:
./ps4_comparesnps.py --nist ${NIST} --snpcalls ${SNPCALLS}
"""

import argparse
import numpy as np
import pandas as pd

def GetCorrect(x):
    if x["nist1"] == x["gt1"] and x["nist2"] == x["gt2"]: return 1
    elif x["nist1"] == x["gt2"] and x["nist2"] == x["gt1"]: return 1
    else: return 0

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--nist", help="NIST calls in tab format", type=str, required=True)
    parser.add_argument("--snpcalls", help="Results of ./ps4_snpcaller.py", type=str, required=True)
    args = parser.parse_args()

    # Load NIST
    nist = pd.read_csv(args.nist, sep="\t")
    nist.columns = ["chrom", "pos", "ref", "nist1", "nist2"]

    # Load snpcalls
    snpcalls = pd.read_csv(args.snpcalls, sep="\t")

    # Merge - only on places where NIST has a call. Use how="outer" to include all positions
    compare = pd.merge(nist, snpcalls, on=["chrom", "pos"], how="left")

    # Overall accuracy rate
    acc = np.mean(compare.apply(GetCorrect, 1))
    print "Overall accuracy at NIST SNP: %s"%acc

    # Accuracy at hets
    acchet = np.mean(compare[compare["nist1"]!=compare["nist2"]].apply(GetCorrect, 1))
    print "Accuracy at hets: %s"%acchet

    # Accuracy at homs
    acchom = np.mean(compare[compare["nist1"]==compare["nist2"]].apply(GetCorrect, 1))
    print "Accuracy at homs: %s"%acchom

    # Acc by coverage
    for covthresh in [5, 10, 20]:
        acccov = np.mean(compare[compare["coverage"]>=covthresh].apply(GetCorrect, 1))
        print "Accuracy at cov >=%s, %s"%(covthresh, acccov)

    # Acc by score
    for sthresh in [1, 2, 30, 50]:
        accscore = np.mean(compare[compare["score"]>=sthresh].apply(GetCorrect, 1))
        print "Accuracy at score >=%s, %s"%(sthresh, accscore)
        
    # Example wrong calls to examine in tview
    correct = compare.apply(GetCorrect, 1)
    print compare[correct==0].head()

if __name__ == "__main__":
    main()
