#!/bin/bash

# Sample annotation and train/test split file (comma-separated text file)
# Rows: samples
# Columns used: sample_id, cancer_status (cancer or normal), train_or_test (train or test)
# This file allows to have other irrelevant columns
input_samples_file="./input/example_multi.cancer.detection_sample.csv"

# Column name and value in "input_samples_file", indicating train/test split.
# For example, "train_or_test:test" indicate Column "train_or_test" for train/test split info and the value in this column "test" means the corresponding row or sample is for testing and the remaining rows/samples for training.
column_test_or_train_info="train_or_test:test"

# Data file  (comma-separated text file)
# Rows: samples
# Columns: sample_id, feature1, feature2, ...
input_data_file="./input/example_multi.cancer.detection_data.csv"

# Output prediction file (comma-separated text file)
out_pred_file="./output/example_multi.cancer.detection_pred.csv"

Rscript ../src/multi.cancer.detection.R \
    $input_samples_file \
    $column_test_or_train_info \
    $input_data_file \
    $out_pred_file
