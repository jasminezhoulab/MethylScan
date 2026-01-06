library(data.table)
library(dplyr)
library(limSolve)
library(RhpcBLASctl)
blas_set_num_threads(2)

source("../src/utility_funs.R")

study_id <- "tissue.deconvolution.study"

args <- commandArgs(trailingOnly = T)

input_markers_file <- args[1]

## Feature matrix file for tissue deconvolution study
## Rows are samples and columns are tissue-specific features
## Columns: sample_id, feature1, feature2, ...
input_feature_matrix_file <- args[2]

## Output tissue deconvolution file name
out_file <- args[3]

## Feature matrix file for tissue deconvolution study
## Rows are samples and columns are tissue-specific features
## Columns: sample_id, feature1, feature2, ...
cat("input_feature_matrix_file\n")
df.data <- read.csv(input_feature_matrix_file, row.names = 1, check.names = FALSE)
df.data[] <- lapply(df.data, function(x) as.numeric(as.character(x)))
df.data <- t(df.data) # transpose to a matrix: rows are features and columns are samples
samples_list <- colnames(df.data)

## Tissue-specific markers file
## Rows are markers, columns are 29 tissues, the value of a marker for a tissue is the reference methylation value of this tissue in this marker.
## Column 1: marker_id
## Column 2+: tissue names
cat("input_markers_file\n")
df.markers <- fread(input_markers_file, header = T, data.table = FALSE)
n_tissues <- ncol(df.markers) - 1
df.markers$marker_id <- as.character(df.markers$marker_id)
rownames(df.markers) <- df.markers$marker_id
df.markers <- df.markers %>% dplyr::filter(marker_id %in% rownames(df.data))
tissues_list <- colnames(df.markers)[2:(1+n_tissues)]

df.data <- df.data[df.markers$marker_id, ]

out_df.tissue.composition <- matrix(NA, nrow = ncol(df.data), ncol = length(tissues_list))
rownames(out_df.tissue.composition) <- samples_list
colnames(out_df.tissue.composition) <- tissues_list

E = rep(1, n_tissues)
F = 1
G = diag(n_tissues)
H = numeric(n_tissues)

A = as.matrix(df.markers[, 2:(1+n_tissues)])

for(sample in samples_list){
  B = df.data[df.markers$marker_id, sample]
  ret <- xsample(A, B, E, F, G, H)
  mean_res <- apply(ret$X, 2, mean)
  out_df.tissue.composition[sample, tissues_list] <- mean_res[tissues_list]
  cat(sample, "\n")
}

tissue_order <- sort(tissues_list)
z<- as.data.frame(cbind(rownames(out_df.tissue.composition), out_df.tissue.composition[, tissue_order]))
colnames(z)[1] <- "sample_id"
fwrite(z, file = out_file, row.names = FALSE, col.names = TRUE, quote=FALSE)
