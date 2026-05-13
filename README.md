
> A high-precision genomic selection pipeline powered by GBLUP, LightGBM, and TabNet. Designed to identify resistant rice accessions using the 44K Rice Diversity Panel (RDP1) and provide ensemble-based resistance scoring.

![R Language](https://img.shields.io/badge/Language-R-276DC3?style=flat-square&logo=r)
![Genomic Selection](https://img.shields.io/badge/Domain-Genomics-success?style=flat-square)
![LightGBM](https://img.shields.io/badge/Model-LightGBM-critical?style=flat-square)
![TabNet](https://img.shields.io/badge/Deep%20Learning-TabNet-blueviolet?style=flat-square)
![GWAS](https://img.shields.io/badge/Analysis-GWAS-orange?style=flat-square)

---

## Overview

The **Rice Blast Resistance Detector (RBRD)** is a comprehensive data science framework built to predict rice blast disease resistance from genomic markers. By integrating traditional quantitative genetics with modern gradient boosting and deep learning, RBRD provides a robust tool for breeders to identify promising resistant varieties without expensive field trials.

The system processes 44,000+ SNP markers across 364 rice accessions, performing automated GWAS to filter significant markers before passing them into a specialized ensemble of GBLUP, LightGBM, and TabNet models.

---

## Features

### Core Pipeline
| Module | Description |
|---|---|
| **GWAS Engine** | Automated Genome-Wide Association Study to identify significant SNP-trait associations. |
| **GBLUP Selection** | Genomic Best Linear Unbiased Prediction using Genomic Relationship Matrices (GRM). |
| **LightGBM Classifier** | Gradient boosting model optimized for large-scale, high-dimensional SNP data. |
| **TabNet Deep Learning** | Attentive tabular deep learning to capture non-linear genetic interactions. |
| **Ensemble Scoring** | Unified **Resistance Score (0-100)** calculated across all three modeling architectures. |

---

## Tech Stack

| Layer | Technology |
|---|---|
| **Programming Language** | R (v4.0+) |
| **Genomic Models** | GBLUP (`rrBLUP`), LightGBM (`lightgbm`), TabNet (`tabnet`) |
| **Statistical Analysis** | GWAS, PCA, Linear Mixed Models |
| **Data Manipulation** | `tidyverse`, `dplyr`, `data.table` |
| **Visualization** | `ggplot2`, `pheatmap`, `pROC` |
| **Machine Learning** | `caret`, `tidymodels`, `torch` (for TabNet) |

---

## Project Structure

```
rice-diversity-project/
├── 1.R                   # Data Loading & Merging (Phenotype + Genotype)
├── 2.R                   # SNP Encoding (Letters to Numbers)
├── 3.R                   # PCA & GWAS Pipeline
├── 4.R                   # Machine Learning Models (LGBM + GBLUP + TabNet)
├── 5.R                   # Final Outputs & Metrics Comparison
├── 6.R                   # Advanced Heatmap Visualizations
├── data/
│   ├── RiceDiversity.44K.MSU6.Genotypes.csv.gz  # DNA Markers
│   ├── RiceDiversity.44K.MSU6.SNP_Information.txt # SNP Annotation
│   └── 12284_2016_116_MOESM1_ESM.xlsx         # Field Trial Scores
└── README.md
```

---

## Getting Started

### 1. Prerequisites
- R (v4.2+ recommended)
- `Rtools` (for Windows users to compile packages)
- Recommended: 16GB+ RAM for handling the 44k SNP matrix

### 2. Setup
```bash
# Clone the repository
git clone https://github.com/Pavan2027/rice-diversity-project.git
cd rice-diversity-project

# Open R or RStudio and install dependencies
Rscript -e "install.packages(c('tidyverse', 'rrBLUP', 'lightgbm', 'tabnet', 'pheatmap', 'pROC', 'caret'))"
```

### 3. Execution Flow
The scripts are designed to be run sequentially:
1. **`1.R`**: Data loading and merging (Phenotype + Genotype).
2. **`2.R`**: SNP encoding and quality control (Letters to Numbers).
3. **`3.R`**: Population structure (PCA) and GWAS pipeline.
4. **`4.R`**: Training Machine Learning models (LGBM, GBLUP, TabNet).
5. **`5.R`**: Performance evaluation and ensemble resistance scoring.
6. **`6.R`**: Advanced heatmap visualizations for genomic similarity.

---

## Results & Performance

The models were evaluated using an 80/20 train-test split, focusing on identifying the "Susceptible" class to ensure no diseased plants are misclassified as resistant.

### **Model Comparison**
| Model | Accuracy | AUC | F1-Score | Kappa |
|---|---|---|---|---|
| **GBLUP** | 0.822 | 0.885 | 0.810 | 0.641 |
| **LightGBM** | **0.863** | **0.912** | **0.855** | **0.724** |
| **TabNet** | 0.849 | 0.898 | 0.841 | 0.698 |

### **Key Insights**
- **Ensemble Robustness**: Combining the models reduces variance, providing a more reliable **Resistance Score (0-100)**.
- **Top Markers**: GWAS identified 19 highly significant SNPs, primarily on Chromosome 6 and 11, which are known hotspots for blast resistance genes (*R-genes*).
- **Population Structure**: PCA reveals that Japonica sub-populations tend to show higher resistance scores in this specific trial.

---

## Implementation Details

### Resistance Scale
The final **Resistance Score** is an ensemble average:
```r
# Resistance Score calculation logic (Simplified)
Resistance_Score = round((1 - (LGBM_Prob + GBLUP_Prob + TabNet_Prob) / 3) * 100, 1)
```
Accessions are classified into five tiers:
- **Highly Resistant** (80-100)
- **Moderately Resistant** (60-80)
- **Intermediate** (40-60)
- **Moderately Susceptible** (20-40)
- **Highly Susceptible** (0-20)

<img width="1350" height="900" alt="resistance_vs_raw_score" src="https://github.com/user-attachments/assets/247da60f-f5ce-4c5d-b7d3-1903de42cc42" />

---

## License

This project is developed for educational purposes (BCSE207L Programming for Data Science). 

**Data Sources:**
- **Phenotypic Data:** [PMC5005242](https://pmc.ncbi.nlm.nih.gov/articles/PMC5005242/) - Genetic Architecture of Rice Resistance to Blast.
- **Genotypic Data:** [Rice Diversity Panel 1 (RDP1)](http://www.ricediversity.org/data/sets/) - 44K SNP Dataset.

---

<p align="center">Built with R + LightGBM + TabNet · 2026</p>
