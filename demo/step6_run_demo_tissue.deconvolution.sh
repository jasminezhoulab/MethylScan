#!/bin/bash

# Marker file with reference methylation levels  (comma-separated text file)
# Rows: markers
# Columns: marker_id, tissue1, tissue2, ...
input_markers_file="./input/example_tissue.deconvolution_marker.reference.matrix.csv"

# Data file  (comma-separated text file)
# Rows: samples
# Columns: sample_id, feature1, feature2, ...
input_sample_data_file="./input/example_tissue.deconvolution_plasma.cfDNA.samples.csv"

# Output tissue fractions file (comma-separated text file)
out_file="./output/example_tissue.deconvolution_tissue.fractions.csv"

Rscript ../src/tissue.deconvolution.study.R \
    $input_markers_file \
    $input_sample_data_file \
    $out_file
