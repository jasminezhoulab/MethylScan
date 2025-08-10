library(data.table)
library(limSolve)
library(RhpcBLASctl)
blas_set_num_threads(2)

source("./common_funs.R")

study_id <- "tissue.deconvolution.study"
n_tissues <- 29

args <- commandArgs(trailingOnly = T)

## Sample annotation file
## Columns: sample_id, disease_type (normal, lung_noncancer_disease, liver_noncancer_disease, lung_cancer, liver_cancer)
input_samples_file <- args[1]

input_markers_file <- args[2]

## Feature matrix file for tissue deconvolution study
## Rows are samples and columns are tissue-specific features
## Columns: sample_id, feature1, feature2, ...
input_feature_matrix_file <- args[2]

## Sample annotation file
## Columns: sample_id, disease_type (normal, lung_noncancer_disease, liver_noncancer_disease, lung_cancer, liver_cancer)
df.samples <- read.csv(input_samples_file)

## Feature matrix file for tissue deconvolution study
## Rows are samples and columns are tissue-specific features
## Columns: sample_id, feature1, feature2, ...
df.data <- read.csv(input_feature_matrix_file, row.names = 1, check.names = FALSE)
df.data[] <- lapply(df.data, function(x) as.numeric(as.character(x)))
df.data <- t(df.data) # transpose to a matrix: rows are features and columns are samples
samples_list <- colnames(df.data)

## Tissue-specific markers file
df.markers <- fread(input_markers_file, header = T, data.table = FALSE)
df.markers$marker_id <- as.character(df.markers$marker_id)
rownames(df.markers) <- df.markers$marker_id
df.markers <- df.markers %>% dplyr::filter(marker_id %in% rownames(df.data))
tissues_list <- colnames(df.markers)[3:2+n_tissues]

df.data <- df.data[df.markers$marker_id, ]

out_df.tissue.composition <- matrix(NA, nrow = ncol(df.data), ncol = length(tissues_list))
rownames(out_df.tissue.composition) <- samples_list
colnames(out_df.tissue.composition) <- tissues_list

E = rep(1, n_tissues)
F = 1
G = diag(n_tissues)
H = numeric(n_tissues)

A = as.matrix(markers.df[, 3:2+n_tissues])

for(sample in samples_list){
  B = all.data[df.markers$marker_id, sample] 
  ret <- xsample(A, B, E, F, G, H)
  mean_res <- apply(ret$X, 2, mean)
  out_df.tissue.composition[sample, tissues_list] <- mean_res[tissues_list]
  cat(sample, "\n")
}

out_dir <- "./output"
dir.create(out_dir, recursive = TRUE)

tissue_order <- sort(tissues_list)
z<- as.data.frame(cbind(rownames(out_df.tissue.composition), out_df.tissue.composition[, tissue_order]))
colnames(z)[1] <- "sample_id"
out_file <- sprintf("%s/%s_tissue.fraction.csv", out_dir, study_id)
fwrite(z, file = out_file, row.names = FALSE, col.names = TRUE, quote=FALSE)
