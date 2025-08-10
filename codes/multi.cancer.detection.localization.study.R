library(dplyr)
library(caret)
library(tidyr)
library(stringr)

source("./common_funs.R")

study_id <- "multi.cancer.detection.localization.study"

four.cancer.types = c(
  "liver_cancer",
  "lung_cancer",
  "ovarian_cancer",
  "stomach_cancer")

args <- commandArgs(trailingOnly = T)

## Sample annocation and 5-fold CV split file
## Columns: sample_id, true_label (cancer or normal), cancer_type (cancer types or normal), fold_id (Fold1, ..., Fold5)
input_samples_cv_split_file <- args[1]

## Feature matrix file for multi-cancer detection
## Rows are samples and columns are pancancer features
## Columns: sample_id, feature1, feature2, ...
input_feature_matrix_file_for_cancer.detection <- args[2]

## Feature matrix file for cancer localization or cancer typing 
## Rows are samples and columns are cancer-typing features
## Columns: sample_id, feature1, feature2, ...
input_feature_matrix_file_for_cancer.typing <- args[3]


classifier_cancer.detection = "LinearSVC_l2_c1"
n.reduced.dim = 500

out_dir <- "./output"
dir.create(out_dir, recursive = TRUE)

## Sample annocation and 5-fold CV split file
## Columns: sample_id, true_label (cancer or normal), cancer_type (cancer types or normal), fold_id (Fold1, ..., Fold5)
df.samples.cv <- read.csv(input_samples_cv_split_file)

## Feature matrix file for multi-cancer detection
## Rows are samples and columns are pancancer features
## Columns: sample_id, feature1, feature2, ...
df.data <- read.csv(input_feature_matrix_file_for_cancer.detection, row.names = 1, check.names = FALSE)
df.data[] <- lapply(df.data, function(x) as.numeric(as.character(x)))

## Perform pan-cancer detection using 5-fold cross validation
fold_ids_list <- sort(unique(df.samples.cv$fold_id))
label_map <- setNames(c(0, 1), c("normal", "cancer"))
all.pred.scores <- c()

for (fold_id in fold_ids_list) {
  train_samples <- df.samples.cv %>% dplyr::filter(fold_id != fold_id)
  train_X <- df.data[train_samples$sample_id, ]
  train_y <- label_map[train_samples$true_label]
  
  test_samples <- df.samples.cv %>% dplyr::filter(fold_id == fold_id)
  test_X <- df.data[test_samples$sample_id, ]
  test_y <- label_map[test_samples$true_label]
  
  z = svd_train_test(train_X, X.test, n.reduced.dim)
  train_X_dim.reduced = z$X_train_reduced
  test_X_dim.reduced = z$X_test_reduced
  
  model <- fit.binary.linear.svm(train_X_dim.reduced, train_y, classifier_cancer.detection)
  pred.scores <- predict.by.binary.linear.svm(model, test_X_dim.reduced)
  names(pred.scores) <- test_samples$sample_id
  all.pred.scores <- c(all.pred.scores, pred.scores)
}

df.samples.cv[names(all.pred.scores), "pred_score_pancancer"] = unname(all.pred.scores)

auc_overall = calculate_auc(label_map[df.samples.cv$true_label],
                                 df.samples.cv$pred_score_pancancer)
result.cancer.detection = calc.auc.sens.spec_cancer.types_and_generate_report_for_cancer.detection(
  df.samples.cv, "cancer_type", "pred_score_pancancer",
  conf.level=NA, max_fp_for_specificity=30, direction="auto")

out_performance_file <- sprintf("%s/multi.cancer.detection.typing.study_performance.csv", out_dir)
cat(sprintf("=== Cancer detection performance ===\nAUC_overall: %.4f\n\n", auc_overall), file = out_performance_file)
write.table(result.cancer.detection, file = out_performance_file,
            sep = ",", row.names = F, col.names = T, append = T)


rm(df.data) # Clear memory

## Perform cancer localization
cancer.detect_fp <- "fp8" # False positive used for obtaining the pan-cancer detection's prediction score cutoff. fp8 corresponds to specificity 98% and means when the number of false positive = 8, we use this corresponding prediction score as the cutoff for binary classification's cutoff.

classifier_cancer.typing = "ovr_LinearSVC_l2_c1"
n.reduced.dim = 300

ratio_cutoff = 1 # ratio of two largest cancer-type membership probabilities, termed as cancer typing confidence score.

## Feature matrix file for cancer localization or cancer typing
## Rows are samples and columns are cancer-typing  features
## Columns: sample_id, feature1, feature2, ...
df.data <- read.csv(input_feature_matrix_file_for_cancer.typing, row.names = 1, check.names = FALSE)
df.data[] <- lapply(df.data, function(x) as.numeric(as.character(x)))

all.pred.scores <- data.frame(matrix(nrow = 0, ncol = 3 + length(four.cancer.types)))
colnames(all.pred.scores) <- c(four.cancer.types, "ratio", "pred.label", "sample_id")

for (fold_id in fold_ids_list) {
  train_samples <- df.samples.cv %>% dplyr::filter(fold_id != fold_id & cancer_type != "normal")
  train_X <- df.data[train_samples$sample_id, ]
  train_y <- train_samples$cancer_type
  
  test_samples <- df.samples.cv %>% dplyr::filter(fold_id == fold_id & cancer_type != "normal")
  test_X <- df.data[test_samples$sample_id, ]
  test_y <- test_samples$cancer_type
  
  z = svd_train_test(train_X, X.test, n.reduced.dim)
  train_X_dim.reduced = z$X_train_reduced
  test_X_dim.reduced = z$X_test_reduced
  
  model <- fit.ovr.linear.svm(train_X_dim.reduced, train_y, classifier_cancer.typing)
  preds <- predict.by.ovr.linear.svm(model, test_X)
  preds$pred.score <- preds$pred.score %>% dplyr::select(all_of(four.cancer.types))
  
  # Generate two additional columns: ratio of two largest prediction scores and pred.label
  cols_pred.scores <- four.cancer.types
  z <- preds$pred.score %>%
    rownames_to_column("sample_id") %>%  # Convert rownames to a column
    rowwise() %>%
    mutate(across(all_of(cols_pred.scores), ~ sigmoid(.))) %>%  # Apply sigmoid
    mutate(
      total = sum(c_across(all_of(cols_pred.scores))),  # Compute sum of transformed scores
      across(all_of(cols_pred.scores), ~ . / total),    # Normalize scores so they sum to 1
      ratio = max(c_across(all_of(cols_pred.scores))) / sort(c_across(all_of(cols_pred.scores)), decreasing = TRUE)[2]  # Compute ratio
    ) %>%
    select(-total) %>%  # Remove the temporary total column
    ungroup() %>%
    column_to_rownames("sample_id") %>%  # Restore row names
    mutate(pred.label = colnames(.)[apply(.[, cols_pred.scores], 1, which.max)]) %>%  # Identify max column name
    mutate(sample_id = rownames(.))
  
  all.pred.scores <- rbind(all.pred.scores, z)
}

## Add cancer typing prediction columns: prediction scores, ratio and pred.label to "df.samples.cv", which are right-adjacent to cancer_type
cols_new <- c(four.cancer.types, "ratio", "pred.label")
df.samples.cv <- df.samples.cv %>%
  left_join(all.pred.scores, by = "sample_id")

# Output prediction scores and true labels for both cancer detection and typing
write.table(df.samples.cv,
            file = sprintf("%s/multi.cancer.detection.typing.study_pred.csv", out_dir),
            sep = ",",
            row.names = F,
            col.names = T)

# Calculate accuracy of cancer typing
cancer.detect_pred.cutoff <- as.numeric(result.cancer.detect[cancer.detect_fp, "threshold"])
cancer.detect_specificity <- result.cancer.detect[cancer.detect_fp, "specificity"]

cat(sprintf("Specificity=%s; %s\n", cancer.detect_fp, cancer.detect_specificity ))

## Get all cancer patients whose cancer.detection prediction score > "cancer.detect_pred.cutoff" (obtained at FP=cancer.detect_fp or specificity=cancer.detect_specificity)
z <- df.samples.cv %>% filter(!!sym(cancer.detect.score_col) > cancer.detect_pred.cutoff)

out_str <- c("",
              sprintf("#All_patients: %d", nrow(df.samples.cv)),
              sprintf("  - #normal (true): %d", nrow(df.samples.cv %>% filter(cancer_type == "normal"))),
              sprintf("  - #cancer (true): %d", nrow(df.samples.cv %>% filter(cancer_type != "normal"))),
              sprintf("#Cancer_patients (true): %d", nrow(df.samples.cv %>% filter(cancer_type != "normal"))), 
              sprintf("#Cancer_patients (pred): %d (pred > %g), at cancer detection's specificity=%s (%s).",
                      nrow(z), cancer.detect_pred.cutoff, cancer.detect_specificity, cancer.detect_fp),
              sprintf("  - #normal (true): %d", nrow(z %>% filter(cancer_type == "normal"))),
              sprintf("  - #cancer (true): %d", nrow(z %>% filter(cancer_type != "normal"))),
              sprintf("#cancer.detect.pred.score_cutoff: %g", cancer.detect_pred.cutoff),
              "")

zz <- z %>% dplyr::filter(ratio >= ratio_cutoff) %>% dplyr::filter(cancer_type != "normal")

ret <- calc_confusion_matrix(zz, "cancer_type", "pred.label", four.cancer.types)

write("\n=== Cancer typing performance ===\n", file = out_performance_file, append = T)
writeLines(out_str,
           file(out_performance_file, open = "a"))
write(ret$csv_lines, file = out_performance_file, append = T)


