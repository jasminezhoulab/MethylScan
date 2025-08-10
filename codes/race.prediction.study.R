library(dplyr)
library(caret)
library(tidyr)
library(stringr)

source("./common_funs.R")

study_id <- "race.prediction.study"

races_list <- c("white", "asian")

args <- commandArgs(trailingOnly = T)

## Sample annocation and 5-fold CV split file
## Columns: sample_id, race, fold_id (Fold1, ..., Fold5)
input_samples_cv_split_file <- args[1]

## Feature matrix file for race prediction study
## Rows are samples and columns are race features
## Columns: sample_id, feature1, feature2, ...
input_feature_matrix_file <- args[2]

classifier = "ovr_LinearSVC_l2_c1"

ratio_cutoff = 1.5 # ratio of two largest race membership probabilities, termed as race prediction confidence score.

## Sample annocation and 5-fold CV split file
## Columns: sample_id, true_label (cancer or normal), cancer_type (cancer types or normal), fold_id (Fold1, ..., Fold5)
df.samples.cv <- read.csv(input_samples_cv_split_file)

## Feature matrix file for race prediction study
## Rows are samples and columns are cancer-typing  features
## Columns: sample_id, feature1, feature2, ...
df.data <- read.csv(input_feature_matrix_file, row.names = 1, check.names = FALSE)
df.data[] <- lapply(df.data, function(x) as.numeric(as.character(x)))

out_dir <- "./output"
dir.create(out_dir, recursive = TRUE)

all.pred.scores <- data.frame(matrix(nrow = 0, ncol = 3 + length(races_list)))
colnames(all.pred.scores) <- c(races_list, "ratio", "pred.label", "sample_id")

for (fold_id in fold_ids_list) {
  train_samples <- df.samples.cv %>% dplyr::filter(fold_id != fold_id)
  train_X <- df.data[train_samples$sample_id, ]
  train_y <- train_samples$race
  
  test_samples <- df.samples.cv %>% dplyr::filter(fold_id == fold_id)
  test_X <- df.data[test_samples$sample_id, ]
  test_y <- test_samples$race
  
  z = svd_train_test(train_X, X.test, n.reduced.dim)
  train_X_dim.reduced = z$X_train_reduced
  test_X_dim.reduced = z$X_test_reduced
  
  model <- fit.ovr.linear.svm(train_X_dim.reduced, train_y, classifier)
  preds <- predict.by.ovr.linear.svm(model, test_X)
  preds$pred.score <- preds$pred.score %>% dplyr::select(all_of(races_list))
  
  # Generate two additional columns: ratio of two largest prediction scores and pred.label
  cols_pred.scores <- races_list
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

## Add cancer typing prediction columns: prediction scores, ratio and pred.label to "df.samples.cv"
cols_new <- c(races_list, "ratio", "pred.label")
df.samples.cv <- df.samples.cv %>%
  left_join(all.pred.scores, by = "sample_id")

# Output prediction scores, predicted labels, and true labels
write.table(df.samples.cv,
            file = sprintf("%s/%s_pred.csv", out_dir, study_id),
            sep = ",", row.names = F, col.names = T)


z <- df.samples.cv %>% dplyr::filter(ratio >= ratio_cutoff)

ret <- calc_confusion_matrix(z, "race", "pred.label", races_list)

out_performance_file <- sprintf("%s/%s_performance.csv", out_dir, study_id)
write("\n=== race prediction performance ===\n", file = out_performance_file, append = T)
writeLines(out_str,
           file(out_performance_file, open = "a"))
write(ret$csv_lines, file = out_performance_file, append = T)
