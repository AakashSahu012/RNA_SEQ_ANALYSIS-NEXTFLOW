# RNA-Seq Analysis Pipeline using Nextflow

This repository contains a reproducible RNA-Seq analysis workflow built with [Nextflow](https://www.nextflow.io/). The pipeline performs quality control, adapter trimming, alignment, quantification, and differential gene expression (DGE) analysis.

## 📁 Workflow Overview

The pipeline consists of the following steps:

1. **Data Download** *(optional)* – SRA Toolkit to download FASTQ files from NCBI.
2. **Quality Control** – FastQC
3. **Trimming** – Trimmomatic
4. **Alignment** – HISAT2
5. **Quantification** – featureCounts
6. **Differential Expression Analysis** – DESeq2
7. **Visualization** – PCA, MA plot, Violin plot

## 🛠️ Tools Used

| Step              | Tool             |
|-------------------|------------------|
| QC                | FastQC           |
| Trimming          | Trimmomatic      |
| Alignment         | HISAT2           |
| Counting          | featureCounts     |
| DGE Analysis      | DESeq2 (in R)    |
| Workflow Mgmt     | Nextflow, Docker |

## 🚀 Getting Started

### Prerequisites

- Nextflow
- Docker or Singularity (for reproducibility)
- (Optional) Conda for R dependencies

### Run the Pipeline

```bash
nextflow run main.nf -c nextflow.config
