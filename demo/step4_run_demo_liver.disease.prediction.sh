#!/bin/bash

# Sample annotation and train/test split file (comma-separated text file)
# Rows: samples
# Columns used: sample_id, disease_type (disease types), train_or_test (train or test)
# This file allows to have other irrelevant columns
input_samples_file="./input/example_liver.disease.prediction_sample.csv"

# Column name and value in "input_samples_file", indicating train/test split.
# For example, "train_or_test:test" indicate Column "train_or_test" for train/test split info and the value in this column "test" means the corresponding row or sample is for testing and the remaining rows/samples for training.
column_test_or_train_info="train_or_test:test"

# Data file  (comma-separated text file)
# Rows: samples
# Columns: sample_id, feature1, feature2, ...
input_data_file="./input/example_liver.disease.prediction_data.csv"

# A list of features selected for liver disease prediction. Each line is a feature name.
input_feature_list_file="./input/example_liver.disease.prediction_feature.txt"

# Output prediction file (comma-separated text file)
out_pred_file="./output/example_liver.disease.prediction_pred.csv"

export PYTHONPATH="$PYTHONPATH:../src"

python ../src/liver.disease.pred.py \
    --sample_split_file=$input_samples_file \
    --sample_column="sample_id" \
    --class_column="disease_type" \
    --column_test_or_train_info=$column_test_or_train_info \
    --input_data_file=$input_data_file \
    --feature_file=$input_feature_list_file \
    --out_pred_file=$out_pred_file
