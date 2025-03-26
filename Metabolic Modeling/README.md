# Li_2025

Metabolic Modeling
Code for metabolic modeling analyses in B Li, 2025

Jason H. Yang Lab @ Rutgers New Jersey Medical School

Author: Gautam Mereddy

These scripts perform metabolic modeling analyses on RNA sequencing profiles generated in B Li, 2025.

Modeling analyses consisted of two steps:
   1. Creating condition-specific models
   2. Creating flux samples

CREATING CONDITION-SPECIFIC MODELS
E coli metabolism was modeled using the iML1515 genome-scale metabolic model (Monk JM, Nat Biotechnol 2017).

RNA sequencing expression profiles were averaged for each biological condition (pEmpty, pF1, pNOX):
- qs-sum-empty.txt
- qs-sum-f1.txt
- qs-sum-nox.txt

Genes were assigned into high (1), medium (0), or low (-1) expression for each biological condition:
- gd-empty.txt
- gd-f1.txt
- gd-nox.txt

iMAT was applied to generate condition-specific models
- iterate_imat.R
- empty.mat
- f1.mat
- nox.mat

FLUX SAMPLING
Flux variability analysis was performed in CobraPy by using optGpSampler. 10,000 flux samples were collected per condition.
- iterate_pysampling.py
- flux-samples-summary.xlsx
