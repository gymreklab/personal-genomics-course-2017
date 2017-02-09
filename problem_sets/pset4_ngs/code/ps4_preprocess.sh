#!/bin/bash

tabix --print-header /oasis/projects/nsf/csd524/mgymrek/data/ps4/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz chr22 | \
    vcf-to-tab | \
    sed 's/INTEGRATION/NIST1\tNIST2/' | sed 's/\//\t/' > \
    /oasis/projects/nsf/csd524/mgymrek/data/ps4/NIST_NA12878_hg38_chr22.tab