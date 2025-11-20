library(dplyr)

source("../../src/utility_funs.R")

load_feature_matrix <- function(file, transpose=TRUE, remove.columns.with.half.NA=TRUE) {
  ## Feature matrix file
  ## Rows are samples and columns are features
  ## Columns: sample_id, feature1, feature2, ...
  if (grepl("\\.gz$", file, ignore.case = TRUE)) {
    df.data <- read.csv(gzfile(file), row.names = 1, check.names = FALSE)
  } else {
    df.data <- read.csv(file, row.names = 1, check.names = FALSE)
  }
  df.data[] <- lapply(df.data, function(x) as.numeric(as.character(x)))
  if (remove.columns.with.half.NA) {
    df.data <- transform_data_cols(df.data)
  }
  if (transpose) {
    return(t(df.data))
  } else {
    return(df.data)
  }
}

args <- commandArgs(trailingOnly = T)

## study_id:
##   - tissue-data-based marker discovery: multi.cancer.study_detection, multi.cancer.study_typing and liver.cancer.study
##   - plasma-data-based marker discovery: liver.disease.study, race.study
study_id <- args[1]
input_samples_file <- args[2]
input_data_matrix_file <- args[3]
out_marker_file <- args[4]

cat(c(sprintf("study_id: %s", study_id),
      sprintf("input_samples_file: %s", input_samples_file),
      sprintf("input_data_matrix_file: %s", input_data_matrix_file),
      sprintf("out_marker_file: %s", out_marker_file)),
      sep = "\n")

if (study_id %in% c("multi.cancer.study_detection",
                    "liver.cancer.study")) {
  diff_cutoff <- as.numeric(args[5])
  topK <- as.numeric(args[6])

  cat(c(sprintf("diff_cutoff: %g", diff_cutoff),
        sprintf("topK_marker: %g", topK)),
        sep = "\n")

  cancer.types_list = c(
    "liver_cancer",
    "lung_cancer",
    "ovarian_cancer",
    "stomach_cancer")
  if (study_id == "liver.cancer.study") {
    cancer.types_list = c("liver_cancer")
  }
  
  ## input_samples_file: rows are tumor samples and columns are sample_id and cancer_type (liver_cancer, lung_cancer, ovarian_cancer, stomach_cancer)
  df_tumor_tissues <- read.csv(input_samples_file, header = T)

  cat(sprintf("#tissue samples for marker discovery: %d\n", nrow(df_tumor_tissues)))
  
  ## input_data_matrix_file: rows are tissue samples and columns are candidate markers
  ## tissue_matrix: rows are candidate markers and columns are samples
  tissue_matrix <- load_feature_matrix(input_data_matrix_file)

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
  writeLines(all.markers, out_marker_file)
  cat(c(sprintf("#cancer_markers: %d", length(all.markers)),
        sprintf("Output: %s", out_marker_file)),
      sep = "\n")
}

if (study_id %in% c("liver.disease.study")) {
  # Column of the test or train set info. For exmaple, an input string "train_or_test:test" indicates all samples (rows) with value "test" of column name "train_or_test" in "input_samples_split_file" are testing samples; and the rest are training samples
  column_test.or.train_info <- args[5] 
  parts <- strsplit(column_test.or.train_info, ":")[[1]]
  column_test.or.train <- parts[1]
  test_value <- parts[2]

  diff_cutoff_for_ovr <- as.numeric(args[6])
  topK_for_ovr <- as.numeric(args[7])

  cat(c(sprintf("diff_cutoff: %g", diff_cutoff_for_ovr),
        sprintf("topK_marker: %g", topK_for_ovr)),
        sep = "\n")

  ## input_samples_file: rows are plasma samples and columns are sample_id and disease_type (liver disease)
  ## Only training samples are used for marker discovery
  train_samples <- read.csv(input_samples_file, header = TRUE) %>%
        dplyr::filter(!!sym(column_test.or.train) != test_value) %>% # Use only training samples
        dplyr::rename(sample = sample_id, type = disease_type)

  cat(sprintf("#train plasma samples for marker discovery: %d\n", nrow(train_samples)))
  
  ## input_data_matrix_file: rows are plasma samples and columns are candidate markers
  ## plasma_matrix: rows are candidate markers and columns are samples
  plasma_matrix <- load_feature_matrix(input_data_matrix_file)
  markers <- select_ovr_features(train_samples, plasma_matrix[, train_samples$sample], cutoff=diff_cutoff_for_ovr, topk.features=topK_for_ovr)
  writeLines(markers, out_marker_file)
  cat(c(sprintf("#liver_disease_markers: %d", length(markers)),
        sprintf("Output: %s", out_marker_file)),
      sep = "\n")
}

if (study_id %in% c("race.study")) {
  # Column of the test or train set info. For exmaple, an input string "train_or_test:test" indicates all samples (rows) with value "test" of column name "train_or_test" in "input_samples_split_file" are testing samples; and the rest are training samples
  column_test.or.train_info <- args[5] 
  parts <- strsplit(column_test.or.train_info, ":")[[1]]
  column_test.or.train <- parts[1]
  test_value <- parts[2]

  diff_cutoff_for_ovr <- as.numeric(args[6])
  topK_for_ovr <- as.numeric(args[7])
  cancer_markers_file <- args[8]

  cancer_markers <- readLines(cancer_markers_file)

  cat(c(sprintf("diff_cutoff: %g", diff_cutoff_for_ovr),
        sprintf("topK_marker: %g", topK_for_ovr),
        sprintf("cancer_markers_file (#marker=%d): %s", length(cancer_markers), cancer_markers_file)),
        sep = "\n")

  ## input_samples_file: rows are plasma samples and columns are sample_id and race
  ## Only training samples are used for marker discovery
  train_samples <- read.csv(input_samples_file, header = TRUE) %>%
        dplyr::filter(!!sym(column_test.or.train) != test_value) %>% # Use only training samples
        dplyr::rename(sample = sample_id, type = race)
  
  cat(sprintf("#train plasma samples for marker discovery: %d\n", nrow(train_samples)))
  
  ## input_data_matrix_file: rows are plasma samples and columns are candidate markers
  plasma_matrix <- load_feature_matrix(input_data_matrix_file)
  ret <- select_ovr_features(train_samples, plasma_matrix[, train_samples$sample], cutoff=diff_cutoff_for_ovr, topk.features=topK_for_ovr)

  # Remove cancer markers from race markers
  markers <- setdiff(ret, cancer_markers)
  writeLines(markers, out_marker_file)
  cat(c(sprintf("#race_markers: %d", length(markers)),
        sprintf("Output: %s", out_marker_file)),
      sep = "\n")
}
