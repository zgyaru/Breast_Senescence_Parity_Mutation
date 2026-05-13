# Breast_Senescence_Parity_Mutation

Time-resolved single-cell atlas of precancerous changes in HR-deficient mouse mammary tissue (Brca1, Brca2, Palb2) across age and parity. Reveals senescence, aberrant differentiation, and immune exhaustion as early hallmarks of breast tumour initiation.

## Overview

This repository contains the analysis code for our study of early tumour-initiating events in homologous recombination (HR)-deficient breast cancer models. We generated a single-cell transcriptomic atlas of >500,000 cells from 62 mice spanning four HR-deficient models (two *Brca1* models, *Brca2*, *Palb2*) and wild-type controls, across a broad age range (up to 2 years) and parity states. Complementary Visium HD spatial transcriptomics was used to map senescent epithelial niches and their immune microenvironment.

Key findings:
- Aberrant differentiation, accelerated ageing, and unresolved senescence of Luminal Adaptive Secretory Precursors (LASPs) across all four HR-deficient models
- Senescent epithelial regions harbour elevated predicted copy number variations and are surrounded by exhausted immune cells
- Pharmacological disruption of senescent cells attenuates aberrant differentiation and reduces the immune-exhausted signature

## Repository structure

| Directory | Contents |
|-----------|----------|
| `preprocessing/` | Raw count QC, doublet removal, ambient RNA correction, and initial filtering
| `scvi/` | scVI integration across mice, genotypes, ages, and parity states
| `sc_annotation/` | Cell type annotation (epithelial, immune, stromal compartments) and marker analysis
| `milo/` | Differential abundance testing across conditions using Milo
| `senescence/` | Senescence scoring, ATB737 treatment analysis and ageing signatures
| `immune_dotplot/` | Immune compartment characterisation and exhaustion signatures
| `Visium_HD/` | Spatial transcriptomics processing, senescent niche mapping, and CNV inference

*Please update the figure mapping above to match your manuscript.*

## Reproducing the analyses

Each subdirectory is self-contained and includes the scripts/notebooks used to generate the corresponding results. 

### Requirements

Analyses were performed primarily in Python (Scanpy, scVI-tools, Milo) and R (inferCNV). Key dependencies:

- Scanpy (v1.11.2);
- anndata (v0.12.10);
- scVI (v1.3.1);
- Scrublet (v0.2.3);
- scipy (v1.16.3);
- Milo (v1.3.1);
- LIANA (v1.5.0);
- LIANA+ (v1.5.0);
- TACCO (v0.4.0);
- infercnv (v3.21);
- Muspan (v1.2.1).

## Data availability

Raw sequencing data will be deposited in a public repository: E-MTAB-16973.
