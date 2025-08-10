# MethylScan codes

## Overview

This code is the implementation of five tasks: (1) multi-cancer detection and localization in the general population, (2) liver cancer surveillance in high-risk individuals, (3) liver disease classification, (4) race prediction, and (5) tissue deconvolution for identification of organ abnormalities.

The first four tasks use Linear Support Vector Machine (LSVM) classifier with the L2 penalty and default hyperparameter C=1 for binary classification (Tasks 1) or multi-class classification. The last task (tissue deconvolution) used the optimization process using R package `limSolve`.

## Prerequisite Packages

This code requires the following R packages:

1) dplyr
2) caret
3) rsample
4) limSolve
5) RhpcBLASctl

## Usage

### Inputs

1) input file of samples' clinical annotation file for multi-cancer detection and typing (including 5-fold cross validation's split), where Columns: sample_id, true_label (cancer or normal), cancer_type (cancer types or normal), fold_id (Fold1, ..., Fold5)
2) input feature matrix file for multi-cancer detection, where Rows are samples and columns are pancancer features and Columns: sample_id, feature1, feature2, ...
3) input feature matrix file for multi-cancer typing, where Rows are samples and columns are cancer-typing features and Columns: sample_id, feature1, feature2, ...
4) input file of samples' clinical annotation file for liver cancer surveillance study (including 5-fold cross validation's split), where sample_id, true_label (cancer or normal), cancer_type (cancer types or normal), fold_id (Fold1, ..., Fold5)
5) input file of samples' clinical annotation file for liver-cancer surveillance (including 5-fold cross validation's split), where sample_id, true_label (cancer or normal), fold_id (Fold1, ..., Fold5)
6) input feature matrix file for liver-cancer surveillance, where Rows are samples and columns are liver-cancer features and Columns: sample_id, feature1, feature2, ...
7) input file of samples' clinical annotation file for liver-etiology prediction (including 5-fold cross validation's split), where columns are sample_id, liver_disease, fold_id (Fold1, ..., Fold5)
8) input feature matrix file for liver-etiology prediction, where Rows are samples and columns are liver-etiology features and Columns: sample_id, feature1, feature2, ...
9) input file of samples' clinical annotation file for race prediction (including 5-fold cross validation's split), where columns are sample_id, race, fold_id (Fold1, ..., Fold5)
10) input feature matrix file for race prediction, where Rows are samples and columns are race-specific features and Columns: sample_id, feature1, feature2, ...
11) input file of samples' clinical annotation file for tissue deconvolution, where columns are sample_id, disease_type (normal, lung_noncancer_disease, liver_noncancer_disease, lung_cancer, liver_cancer)
12) input feature matrix file for tissue deconvolution, where Rows are samples and columns are tissue-specific features and Columns: sample_id, feature1, feature2, ...
13) input tissue-specific markers file for tissue deconvolution, where Rows are markers, columns are 29 tissues, the value of a marker for a tissue is the reference methylation value of this tissue in this marker.
