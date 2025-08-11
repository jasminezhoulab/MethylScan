library(dplyr)

source("./common_funs.R")

args <- commandArgs(trailingOnly = T)

## study_id:
##   - tissue-data-based marker discovery: multi.cancer.detection.localization.study and liver.cancer.surveillance.study
##   - plasma-data-based marker discovery: liver.etiology.prediction.study, race.prediction.study
study_id <- args[1]
input_samples_file <- args[2]
input_data_matrix_file <- args[3]

out_dir <- "./output"
dir.create(out_dir, recursive = T)

if (study_id %in% c("multi.cancer.detection.localization.study",
                    "liver.cancer.surveillance.study")) {
  diff_cutoff <- as.numeric(args[4])
  topK <- as.numeric(args[5])
  cancer.types_list = c(
    "liver_cancer",
    "lung_cancer",
    "ovarian_cancer",
    "stomach_cancer")
  if (study_id == "liver.cancer.surveillance.study") {
    cancer.types_list = c("liver_cancer")
  }
  
  ## input_samples_file: rows are tumor samples and columns are sample_id and cancer_type (liver_cancer, lung_cancer, ovarian_cancer, stomach_cancer)
  df_tumor_tissues <- read.csv(input_samples_file, header = T)
  
  ## input_data_matrix_file: rows are candidate markers and columns are all tissue samples
  if (grepl("\\.gz$", input_data_matrix_file)) {
    tissue_matrix <- read.csv(gzfile(input_data_matrix_file), check.names = FALSE, row.names = 1)
  } else {
    tissue_matrix <- read.csv(input_data_matrix_file, check.names = FALSE, row.names = 1)
  }
  
  all.markers <- c()
  for (cancer_type in cancer.types_list) {
    tumor_tissues_list <- df_tumor_tissues %>% dplyr::filter(cancer_type == cancer_type)
    markers <- select_markers_by_paired.tissues_and_occurrence(
      tumor_tissues_list$sample_id,
      tissue_matrix,
      diff_cutoff,
      topK)
    all.markers <- c(all.markers, markers)
  }
  all.markers <- unique(all.markers)
  
  writeLines(all.markers,
             sprintf("%s/%s_markers.txt", out_dir, study_id))
}

if (study_id %in% c("liver.etiology.prediction.study")) {
  diff_cutoff_for_ovr <- as.numeric(args[4])
  topK_for_ovr <- as.numeric(args[5])
  diff_cutoff_for_ovo <- as.numeric(args[6])
  topK_for_ovo <- as.numeric(args[7])
  ## input_samples_file: rows are plasma samples and columns are sample_id and class_names (either liver diseases or races)
  df_plasma_samples <- read.csv(input_samples_file, header = T)
  df_plasma_samples <- df_plasma_samples %>% dplyr::rename(sample = sample_id, type = liver_disease)
  
  ## input_data_matrix_file: rows are candidate markers and columns are all plasma samples
  if (grepl("\\.gz$", input_data_matrix_file)) {
    plama_matrix <- read.csv(gzfile(input_data_matrix_file), check.names = FALSE, row.names = 1)
  } else {
    plama_matrix <- read.csv(input_data_matrix_file, check.names = FALSE, row.names = 1)
  }
  
  ret1 <- select_ovr_features_by.limma(df_plasma_samples, plama_matrix, cutoff=diff_cutoff_for_ovr, topk.features=topK_for_ovr)
  ret2 <- select_ovr_features_by.limma(df_plasma_samples, plama_matrix, cutoff=diff_cutoff_for_ovo, topk.features=topK_for_ovo)
  markers <- unique(c(ret1$union, ret2$union))
  return(markers)
}

if (study_id %in% c("race.prediction.study")) {
  diff_cutoff_for_ovr <- as.numeric(args[4])
  topK_for_ovr <- as.numeric(args[5])
  ## input_samples_file: rows are plasma samples and columns are sample_id and class_names (either liver diseases or races)
  df_plasma_samples <- read.csv(input_samples_file, header = T)
  df_plasma_samples <- df_plasma_samples %>% dplyr::rename(sample = sample_id, type = race)
  
  ## input_data_matrix_file: rows are candidate markers and columns are all plasma samples
  if (grepl("\\.gz$", input_data_matrix_file)) {
    plama_matrix <- read.csv(gzfile(input_data_matrix_file), check.names = FALSE, row.names = 1)
  } else {
    plama_matrix <- read.csv(input_data_matrix_file, check.names = FALSE, row.names = 1)
  }
  
  ret <- select_ovr_features_by.limma(df_plasma_samples, plama_matrix, cutoff=diff_cutoff_for_ovr, topk.features=topK_for_ovr)
  markers <- ret$union
  return(markers)
}
