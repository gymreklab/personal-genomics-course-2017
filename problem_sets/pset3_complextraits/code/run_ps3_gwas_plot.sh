#!/bin/bash

USERPREFIX=/oasis/projects/nsf/csd524/$USER/ps3/results/ps3_gwas

# Manhattan plot
./ps3_manhattan.py \
    --assoc ${USERPREFIX}.assoc.linear.tab \
    --out ${USERPREFIX}.manhattan.png
./ps3_manhattan.py \
    --assoc ${USERPREFIX}.withcovar.assoc.linear.tab \
    --out ${USERPREFIX}.withcovar.manhattan.png

# QQ plot
./ps3_qqplot.py \
    --assoc ${USERPREFIX}.assoc.linear.tab \
    --out ${USERPREFIX}.qq.png
./ps3_qqplot.py \
    --assoc ${USERPREFIX}.withcovar.assoc.linear.tab \
    --out ${USERPREFIX}.withcovar.qq.png

# PCA
./ps3_plotpcs.py \
    --pcs ${USERPREFIX}.eigenvec \
    --out ${USERPREFIX}.pcplot.png