# MethylScan Codes

## Overview

This repository contains the code, pipeline, and demo for the five studies of **MethylScan**: multi-cancer detection and typing (i.e., tissue-of-origin), liver cancer surveillance (i.e., classifying between liver cancer and high-risk patients), non-cancer liver disease classification, tissue deconvolution, and race prediction.

---
## Prerequisite Packages
All the studies were developed in on UCLA Hoffman2 cluster. All the analyses in this study were done there. The system information of UCLA Hoffman2 cluster is:
Linux n6015 2.6.32-754.14.2.el6.x86_64 #1 SMP Tue May 14 19:35:42 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux. There is no non-standard hardware required.

The following R packages are required:

- `R 4.2.2`
- `dplyr 1.1.4` 
- `LiblineaR 2.10.23` 
- `caret 6.0.94`  
- `rsample 1.2.1`  
- `limSolve 1.5.7.1`  
- `RhpcBLASctl 0.23.42`  
- `limma 3.54.2`
- `matrixStats 1.2.0`
- `tidyr 1.3.1`
- `stringr 1.5.1`
- `datatable 1.15.2`
- `rlang 1.1.3`
- `tibble 3.2.1`

The following Python packages are required:

- `Python 3.8.5`
- `numpy 1.24.4`  
- `scikit-learn 1.4.dev0`  
- `pandas 2.0.3`
- `csv 1.0`

## Installation

All codes were implemented using R and python scripts. There is no specific installation needed, if all required packages and tools listed in PREREQUISITE PACKAGES are successfully installed. All the source codes are saved under src folder.

Although the R and python code was developed in the environment of UCLA Hoffman2 cluster, the codes and demo bash scripts can be directly run in other Linux environment.

---
## Repository Structure

- Codes: `src/`. Contains R and python codes for all studies. Each code file typically:
  - Loads the data matrix and sample annotations (including true labels and training/testing split information)   
  - Trains the predictive model
  - Make predictions and output
- Demo data and pipeline: `demo/`.

---
## Usage

## Input Parameters and File Formats

#### **Multi-Cancer Detection** (inputs of `multi.cancer.detection.R`)
- input sample info file, including columns regarding to clinical annotations, i.e., `sample_id`, `cancer_status` (normal or cancer) and training/testing sample split info, i.e., column_test.or.train (value is a string that indicates test samples/rows)
- input data matrix file:
  - Rows are samples and columns are cancer-specific features  
  - Columns: sample, feature1, feature2, ...  
- a string <column_name:test_sample_indicator_string>, where "column_name" is the column name of indicating a row representing a training or testing sample and "test_sample_indicator_string" represents the value of testing sample's indicator. For example, for the input string "train_or_test:test", it indicates all samples (rows) with value "test" of column name "train_or_test" in the input sample info file are testing samples; and the rest are training samples
- output prediction file, including all columns of the input sample info file and the additional column "pred.score_multi.cancer".
- Note:
  - The classifier Linear SVM with default parameters and the dimension reduction have been embedded in the code.

#### **Multi-Cancer Typing** (inputs of `multi.cancer.typing.R`)
- input sample info file, including columns regarding to clinical annotations, i.e., `sample_id`, `cancer_status` (normal or cancer), `cancer_type`, and training/testing sample split info, i.e., column_test.or.train (value is a string that indicates test samples/rows)
- input data matrix file:
  - Rows are samples and columns are cancer-type-specific features  
  - Columns: sample, feature1, feature2, ...  
- a string <column_name:test_sample_indicator_string>, where "column_name" is the column name of indicating a row representing a training or testing sample and "test_sample_indicator_string" represents the value of testing sample's indicator. For example, for the input string "train_or_test:test", it indicates all samples (rows) with value "test" of column name "train_or_test" in the input sample info file are testing samples; and the rest are training samples
- output prediction file, including all columns of the input sample info file and the additional columns "cancer_type_1", "cancer_type_2", ..., "ratio", "pred_label", where the columns "cancer_type_*" are the cancer-type-specific probability of each cancer type and their sum equals to one, the column "ratio" is the ratio of the largest and second-largest cancer type probablities, and the column "pred_label" is the predicted cancer type that has the largest probability.
- Note:
  - The multi-class classifier Linear SVM with default parameters and the dimension reduction have been embedded in the code.

#### **Liver Cancer Surveillance** (inputs of `liver.cancer.surveillance.R`)
- input sample info file, including columns regarding to clinical annotations, i.e., `sample_id`, `cancer_status` (normal or cancer) and training/testing sample split info, i.e., column_test.or.train (value is a string that indicates test samples/rows)
- input data matrix file:
  - Rows are samples and columns are liver-cancer-specific features  
  - Columns: sample, feature1, feature2, ...  
- a string <column_name:test_sample_indicator_string>, where "column_name" is the column name of indicating a row representing a training or testing sample and "test_sample_indicator_string" represents the value of testing sample's indicator. For example, for the input string "train_or_test:test", it indicates all samples (rows) with value "test" of column name "train_or_test" in the input sample info file are testing samples; and the rest are training samples
- output prediction file, including all columns of the input sample info file and the additional column "pred.score_liver.cancer".
- Note:
  - The classifier Linear SVM with default parameters and the dimension reduction have been embedded in the code.

#### **Liver Disease Classification** (inputs of `liver.disease.pred.py`)
- input sample info file, including columns regarding to clinical annotations (sample ID column and liver disease type column), and training/testing sample split info, i.e., column_test.or.train (value is a string that indicates test samples/rows)
- a string <sample_column>, indicating sample ID's column name
- a string <class_column>, indicating liver disease type's column name
- input data matrix file:
  - Rows are samples and columns are features (not necessarily liver-disease-specific)
  - Columns: sample, feature1, feature2, ... 
- input liver-disease-specific features file: each line is a feature name. 
- a string <column_name:test_sample_indicator_string>, where "column_name" is the column name of indicating a row representing a training or testing sample and "test_sample_indicator_string" represents the value of testing sample's indicator. For example, for the input string "train_or_test:test", it indicates all samples (rows) with value "test" of column name "train_or_test" in the input sample info file are testing samples; and the rest are training samples
- output prediction file, including all columns of the input sample info file and the additional columns: The column "pred_label" is the predicted liver disease that has the largest probability. The column "ratio" is the ratio of the largest and second-largest liver disease type probablities. The columns "liver_disease_type_1", "liver_disease_type_2", ..., "ratio", "pred_label", where the columns "liver_disease_type_*" are the liver-disease-type-specific probability of each liver disease and their sum equals to one. 
- Note:
  - The multi-class classifier Linear SVM with default parameters have been embedded in the code.

#### **Race Prediction** (inputs of `liver.disease.pred.py`)
- input sample info file, including columns regarding to clinical annotations, i.e., `sample_id`, `race`, and training/testing sample split info, i.e., column_test.or.train (value is a string that indicates test samples/rows)
- a string <column_name:test_sample_indicator_string>, where "column_name" is the column name of indicating a row representing a training or testing sample and "test_sample_indicator_string" represents the value of testing sample's indicator. For example, for the input string "train_or_test:test", it indicates all samples (rows) with value "test" of column name "train_or_test" in the input sample info file are testing samples; and the rest are training samples
- input data matrix file:
  - Rows are samples and columns are features (not necessarily race-specific)
  - Columns: sample, feature1, feature2, ...  
- input race-specific features file: each line is a feature name.
- output prediction file, including all columns of the input sample info file and the additional columns "race_1", "race_2", ..., "ratio", "pred_label", where the columns "race_*" are the race-specific probability of each cancer type and their sum equals to one, the column "ratio" is the ratio of the largest and second-largest race probablities, and the column "pred_label" is the predicted race that has the largest probability.
- Note:
  - The multi-class classifier Linear SVM with default parameters and the dimension reduction have been embedded in the code.


---

## Step-by-step pipeline with demo of example data

The command-line arguments, and detailed usage instructions are provided in the folder `demo/`

### Command line scripts that run demo data for each study

- Step 0: Ensure python and its packages to be installed. Script is `step0_ensure_R_and_libraries.sh`. Note that R and its libraries need to be installed manually.
- Step 1: Run multi-cancer detection study with the demo data. Script is `step1_run_demo_multi.cancer.detection.sh`.
- Step 2: Run multi-cancer typing study with the demo data. Script is `step2_run_demo_multi.cancer.typing.sh`.
- Step 3: Run liver-cancer surveillance study with the demo data. Script is `step3_run_demo_liver.cancer.surveillance.sh`.
- Step 4: Run liver disease prediction study with the demo data. Script is `step4_run_demo_liver.disease.prediction.sh`.
- Step 5: Run race prediction study with the demo data. Script is `step5_run_demo_race.prediction.sh`.
- Step 6: Run tissue deconvolution study with the demo data. Script is `step6_run_demo_tissue.deconvolution.sh`.

### Demo data

A demo is provided along with the scripts. Example data and required reference files are saved under demo folder.

This demo can be run with command lines. Please run the following commands:

```
cd demo
./step0_ensure_R_and_libraries.sh
./step1_run_demo_multi.cancer.detection.sh
./step2_run_demo_multi.cancer.typing.sh
./step3_run_demo_liver.cancer.surveillance.sh
./step4_run_demo_liver.disease.prediction.sh
./step5_run_demo_race.prediction.sh
./step6_run_demo_tissue.deconvolution.sh
```

In this demo, six output prediction files are generated:

1) output/example_multi.cancer.detection_pred.csv
2) output/example_multi.cancer.typing_pred.csv
3) output/example_liver.cancer.surveillance_pred.csv
4) output/example_liver.disease.prediction_pred.csv
5) output/example_race_pred.csv
6) output/example_tissue.deconv_pred.csv

Please compare each output file with its reference file, which is already provided in the package:

1) output/example_multi.cancer.detection_pred.csv.reference
2) output/example_multi.cancer.typing_pred.csv.reference
3) output/example_liver.cancer.surveillance_pred.csv.reference
4) output/example_liver.disease.prediction_pred.csv.reference
5) output/example_race_pred.csv.reference
6) output/example_tissue.deconv_pred.csv.reference
