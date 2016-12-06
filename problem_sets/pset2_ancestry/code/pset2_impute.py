#!/usr/bin/env python

"""
CSE291 - Personal Genomics for Bioinformaticians
Problem Set 2 - Ancestry: Imputations

./pset2_impute.py \
  /oasis/projects/nsf/csd524/$USER/ps2/results/ps2_impute.combined \
  /oasis/projects/nsf/csd524/mgymrek/data/ps2/ps2_impute.heldout.gen.gz \
  /oasis/projects/nsf/csd524/mgymrek/data/ps2/imputation/1000GP_Phase3/1000GP_Phase3_chr16.legend.gz
"""

import pandas as pd
from scipy.stats import pearsonr
import sys

try:
    imputefile = sys.argv[1]
    truthfile = sys.argv[2]
    legendfile = sys.argv[3]
except:
    sys.stderr.write(__doc__+"\n")
    sys.exit(1)

# Load imputation results
gts = ["00","01","11"]
header = ["ident", "position", "ref1", "ref2"] + \
    ["ceu_%s"%gt for gt in gts] + \
    ["yri_%s"%gt for gt in gts] + \
    ["ceu_yri_%s"%gt for gt in gts] + \
    ["not_asw_%s"%gt for gt in gts]
impres = pd.read_csv(imputefile, sep=" ", names=header)

# Load truth
truth = pd.read_csv(truthfile, sep=" ", names=["chrom","rsid","position","ref","alt","truth_00", "truth_01", "truth_11"])
impres = pd.merge(impres, truth, on=["position"])

# Load legend
legend = pd.read_csv(legendfile, sep=" ")

# merge
data = pd.merge(impres, legend, on=["position"])

# Annotate best gt and score for each panel
def GetGenotype(x00, x01, x11):
    genotypes = [0, 1, 2]
    scores = [x00, x01, x11]
    ind = scores.index(max(scores))
    return genotypes[ind]

data["truth_gt"] = data.apply(lambda x: GetGenotype(x["truth_00"], x["truth_01"], x["truth_11"]), 1)
for refpanel in ["ceu", "yri", "ceu_yri", "not_asw"]:
    data["%s_gt"%refpanel] = data.apply(lambda x: GetGenotype(x["%s_00"%refpanel], x["%s_01"%refpanel], x["%s_11"%refpanel]), 1)
    data["%s_score"%refpanel] = data.apply(lambda x: max([x["%s_00"%refpanel], x["%s_01"%refpanel],x["%s_11"%refpanel]]), 1)
    print refpanel, pearsonr(data["truth_gt"], data["%s_gt"%refpanel]), data.shape[0]
