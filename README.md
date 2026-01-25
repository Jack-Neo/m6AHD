
# m6AHD: A Framework for Identifying Abnormal m6A Modifications in Heart Diseases

**m6AHD** (m6A in Heart Disease) is a bioinformatics framework designed to predict aberrant N6-methyladenosine (m6A) modification sites associated with various cardiac pathological states. By integrating multi-condition MeRIP-seq data and sequence-derived features, this tool utilizes Random Forest algorithms to identify dysregulated (upregulated or downregulated) m6A sites.

This repository contains the source code for data preprocessing, feature extraction, and model construction/prediction as described in our manuscript.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Repository Structure](#repository-structure)
- [Usage](#usage)
  - [1. Data Preprocessing](#1-data-preprocessing)
  - [2. Model Prediction](#2-model-prediction)
- [Citation](#citation)


## Prerequisites

To run the scripts in this repository, you need **R** installed (tested on version 4.4.2). 
Please ensure the following R packages are installed:

```r
install.packages(c("randomForest", "caret", "pROC", "Biostrings", "seqinr"))
# Note: List any other specific packages your code uses
```
## Repository Structure
```text
m6AHD/
├── data_site.R               # Script for identifying differentially methylated sites
├── final model/              # Contains the specific models for 5 cardiac conditions
│   ├── AC16EV/               # Evodiamine-induced cardiotoxicity
│   │   ├── up/               # Model for Upregulated sites
│   │   │   └── final model.R
│   │   └── down/             # Model for Downregulated sites
│   │       └── final model.R
│   ├── AC16Mat/              # Matrine-induced cardiotoxicity
│   │   ├── up/
│   │   └── down/
│   ├── Calcific/             # Heart calcification
│   │   ├── up/
│   │   └── down/
│   ├── hypertrophy/          # Cardiac hypertrophy
│   │   ├── up/
│   │   └── down/
│   └── TKIs/                 # TKI-induced cardiotoxicity
│       ├── up/
│       └── down/
│
└──train
│
└── README.md
```
## Usage
**1**. Data Preprocessing
The script data_site.R is the starting point of the pipeline. Its primary function is to classify m6A sites into Upregulated or Downregulated categories based on the comparison between the experimental (disease) group and the control group.

Functionality:

- Inputs: MeRIP-seq peak calling results or methylation level matrices.

- Logic:

  - Upregulated: Sites present/methylated in the experimental group but absent/hypomethylated in the control group.

  - Downregulated: Sites absent/hypomethylated in the experimental group but present/methylated in the control group.
 
 How to run:
```r
source("data_site.R")
# Ensure your input data is placed in the working directory as required by the script.
``` 
**2**. Model Prediction
    The `final model` directory contains the optimized Random Forest models for five specific cardiac pathological conditions. Each condition is further divided into `up` (for predicting hypermethylation) and `down` (for predicting hypomethylation).
**Pathological Conditions:**

1.  **AC16EV:** Evodiamine-induced cardiotoxicity
    
2.  **AC16Mat:** Matrine-induced cardiotoxicity
    
3.  **TKIs:** TKI-induced cardiotoxicity
    
4.  **hypertrophy:** Cardiac hypertrophy
    
5.  **Calcific:** Heart calcification
    

**Running a Specific Model:** Each sub-folder (e.g., `final model/AC16EV/up/`) contains a script named `final model.R`. This script loads the training data, processes sequence features (e.g., One-Hot, NAC, PseKNC), and executes the Random Forest algorithm with optimized parameters (ntree, mtry).

**Example:** To run the model for Evodiamine-induced upregulation:

1.  Navigate to the specific directory: `final model/AC16EV/up/`
    
2.  Run the script:
```r
setwd("./final model/AC16EV/up/")
source("final model.R")
``` 
Note: Ensure your input data is placed in the working directory as required by the script. The `final model.R` script will output the performance metrics (AUROC).



## Citation
If you use this code or data in your research, please cite our paper:

> **m6AHD: A new framework for identifying abnormal N6-Methyladenosine (m6A) in heart diseases based on sequencing features**  Jiajie Lu, Yanan Li, Yuxiang Hong, Dongshan Liao, Guanhua Fang)
