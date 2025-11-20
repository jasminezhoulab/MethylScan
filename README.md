# MethylScan Codes

## Overview

This repository contains the analysis code for **MethylScan™**, a cell-free DNA (cfDNA) methylome–based framework for cancer early detection, disease classification, and tissue-signal interpretation.

MethylScan supports **five analysis tasks**:

1. **Multi-cancer detection and tissue-of-origin localization**  
2. **Liver cancer surveillance** 
3. **Liver disease classification** 
4. **Race prediction**  
5. **Tissue deconvolution**

Tasks 1–4 use **Linear Support Vector Machine (LSVM)** models with L2 regularization and hyperparameter **C = 1**.  
Task 5 uses constrained optimization via **`limSolve`**.

Most code is written in **R**, with optional helper scripts in **Python**.

---

## Step-by-Step Pipeline

We have provided a comprehensive **step-by-step pipeline**, **example datasets**, and **detailed CLI documentation** in the **GitHub link referenced in the manuscript**.  
Please refer to that page for exact replication of the analyses described in the paper.

This README provides a high-level overview of directory structure, input formats, and usage.

---

## Required R Packages

The following R packages are required:

- `dplyr`  
- `caret`  
- `rsample`  
- `limSolve`  
- `RhpcBLASctl`  
- `limma`

Optional R packages:

- `ggplot2`  
- `pROC` / `PRROC`  
- `readr`, `data.table`, `vroom`  

Optional Python (for helper tools):

- `numpy`  
- `scikit-learn`  
- `scipy`

---

## Repository Structure

### `src/` — Core Analysis Code

Contains R scripts for all five MethylScan tasks:

- **Task 1: Multi-cancer detection & localization**  
- **Task 2: Liver cancer surveillance**  
- **Task 3: Liver disease classification**  
- **Task 4: Race prediction**  
- **Task 5: Tissue deconvolution**

Each script typically:

- Loads feature matrices & clinical annotations  
- Applies 5-fold cross-validation using fold IDs  
- Trains models (LSVM or constrained solving)  
- Outputs predictions and performance summaries  

### `demo/` — Example Data & Templates

Contains:

- Example annotation files  
- Example feature matrices  
- Tissue reference marker files (for Task 5)  
- Toy data for demonstration  
- Template input formats for all five tasks

Use these as a reference for formatting your own cohort data.

---

## Input File Formats

All input files are tab-delimited with headers.

### Common Requirements Across Tasks

- `sample_id`  
- `fold_id` (Fold1–Fold5)  
- Task-specific label column  
- Feature matrix with columns:  
  - `sample_id`, `feature1`, `feature2`, …  

### Task-specific Inputs

#### **1. Multi-Cancer Detection & Localization**
- Clinical annotation (`true_label`, `cancer_type`, `fold_id`)
- Feature matrices:
  - Detection features  
  - Typing (tissue-of-origin) features  

#### **2. Liver Cancer Surveillance**
- Annotation with:
  - `true_label`  
  - `cancer_type` or binary label  
  - `fold_id`
- Feature matrix for liver surveillance

#### **3. Liver Disease Classification**
- Annotation: `liver_disease`, `fold_id`
- Feature matrix

#### **4. Race Prediction**
- Annotation: `race`, `fold_id`
- Feature matrix

#### **5. Tissue Deconvolution**
- Annotation: `disease_type`
- cfDNA marker matrix
- 29-tissue reference marker matrix

---

## Usage

The exact command-line arguments, flags, and detailed usage instructions are provided in the **step-by-step pipeline webpage linked in the manuscript**.

Below is a general example:

```bash
Rscript src/run_multi_cancer_detection.R \
  --clinical demo/multi_cancer/clinical_annotation.tsv \
  --feature_detection demo/multi_cancer/features_detection.tsv \
  --feature_typing demo/multi_cancer/features_typing.tsv \
  --out_dir results/multi_cancer/

## Demo
The demo/ directory includes:
- Example sample annotation and split file for each task
- Example data file for each task
