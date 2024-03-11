---
title: SNPAssignment
description: "The script used to create .annot files for the data."
sidebar_position: 2
---

# SNPAssignment

## Explanaiton

In order to use partitioned heritability with LDSC, an annotation file is required. An annotation file is just a classification of every SNP in the given PLINK files. So if a state spans the region 400-800 on chromosome 2, then all SNPs that lie in this genomic region are classified as being in that state. The output `.annot` file that comes out of this R script is effectively a `.bim` file with additional Boolean columns signifying which state each SNP belongs to.

## Overlaps

If one looks at the official documentation for LDSC, they will see information given for  overlapping annotations, frq files and a 'base' characterisation of SNPs. All of these are not actually required in the case of ChromOptimise, so they are ignored. The reason for this is due to the fact that our state assignments (extracted from ChromHMM) necessarily will not overlap. A genomic region can only be assigned exactly one state (never more, never less), as such there will be no overlapping regions and so the extra flags in the ldsc command are therefore not required.

This also has the added benefit of making our `.annot` files having no colinearity between annoation columns (so ldsc should never fail from colinearity errors).