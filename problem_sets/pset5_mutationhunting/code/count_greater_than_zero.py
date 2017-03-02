#!/usr/bin/env python

"""
Utility script to count all columns with a value > 0

Usage:
./count_greater_than_zero.py file startcol endcol
"""

import sys

try:
    infile = sys.argv[1]
    startcol = int(sys.argv[2])
    endcol = int(sys.argv[3])
except:
    sys.stderr.write(__doc___)
    sys.exit(1)

if infile == "-": inf = sys.stdin
else: inf = open(infile, "r")

line = inf.readline()
while line != "":
    # TODO process line want count of > 0, not cols themselves
    outitems = []
    items = line.strip().split()
    g0 = []
    for i in range(len(items)):
        if i < startcol or i > endcol:
            outitems.append(items[i])
        else: g0.append(int(int(items[i])>0))
    if sum(g0) > 0:
        sys.stdout.write("\t".join(map(str, outitems+[sum(g0)]))+"\n")
    line = inf.readline()
