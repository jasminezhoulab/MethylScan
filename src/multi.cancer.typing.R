library(dplyr)
library(rlang)
library(tibble)

source("../src/utility_funs.R")

study_id <- "multi.cancer.study_typing"

args <- commandArgs(trailingOnly = T)

## Sample annotation and train/test split file
## Columns: sample_id, cancer_status (cancer or normal), cancer_type (cancer types or normal), column_test.or.train (value is a string that indicates test samples/rows)
input_samples_split_file <- args[1]

# Column of the test or train set info. For exmaple, an input string "train_or_test:test" indicates all samples (rows) with value "test" of column name "train_or_test" in "input_samples_split_file" are testing samples; and the rest are training samples
column_test.or.train_info <- args[2] 
parts <- strsplit(column_test.or.train_info, ":")[[1]]
column_test.or.train <- parts[1]
test_value <- parts[2]

## Feature matrix file for cancer localization or cancer typing 
## Rows are samples and columns are cancer-typing features
## Columns: sample_id, feature1, feature2, ...
input_feature_matrix_file_for_cancer.typing <- args[3]

## Output prediction file name
out_file <- args[4]

classifier_cancer.typing = "ovr_LinearSVC_l2_c1" # One-vs-Rest Linear SVM with L2 regularization and Cost=1
n.reduced.dim_cancer.typing = 300

## Sample annotations and train/test split file
## Columns: sample_id, cancer_status (cancer or normal), cancer_type (cancer types or normal), column_test.or.train (train or test)
df.samples <- read.csv(input_samples_split_file, header=TRUE)
rownames(df.samples) <- df.samples$sample_id

cancer.types <- sort(unique(df.samples$cancer_type))

cat(c(sprintf("Study: %s (%s)", study_id, paste(cancer.types, collapse = ",")),
      sprintf("Sample split file: %s", input_samples_split_file),
      sprintf("Column (test or train indicator): %s (%s)", column_test.or.train, test_value)),
      sprintf("Feature matrix file (row as sample; column as feature): %s", input_feature_matrix_file_for_cancer.typing),
      sep = "\n")

########################################################################
## Perform cancer typing

## Feature matrix file for cancer localization or cancer typing
## Rows are samples and columns are cancer-typing  features
## Columns: sample_id, feature1, feature2, ...
df.data <- read.csv(input_feature_matrix_file_for_cancer.typing, row.names = 1, check.names = FALSE)
df.data[] <- lapply(df.data, function(x) as.numeric(as.character(x)))

cat(c(sprintf("#markers: %d", ncol(df.data)),
      sprintf("Classifier: %s", classifier_cancer.typing)),
      sep = "\n")

train_samples <- df.samples %>% dplyr::filter(!!sym(column_test.or.train) != test_value & cancer_status != "normal")
train_X <- df.data[train_samples$sample_id, ]
train_y <- train_samples$cancer_type
counts <- sapply(cancer.types, \(x) sum(train_samples$cancer_type == x))
cat(sprintf("#train_set: %d (%s)\n", nrow(train_X), paste(sprintf("%s: %d", names(counts), counts), collapse = ", ")))

test_samples_for_typing <- df.samples %>% dplyr::filter(!!sym(column_test.or.train) == test_value & cancer_status != "normal")
test_X <- df.data[test_samples_for_typing$sample_id, ]
test_y <- test_samples_for_typing$cancer_type
counts <- sapply(cancer.types, \(x) sum(test_samples_for_typing$cancer_type == x))
cat(sprintf("#test_set: %d (%s)\n", nrow(test_X), paste(sprintf("%s: %d", names(counts), counts), collapse = ", ")))

if (n.reduced.dim_cancer.typing >= ncol(train_X)) {
  train_X_dim.reduced = train_X
  test_X_dim.reduced = test_X
  cat(sprintf("DO NOT perform dimension reduction, because #marker=%d <= reduced_dimension=%d\n", ncol(train_X), n.reduced.dim_cancer.typing))
} else {
  z = svd_train_test(train_X, test_X, n.reduced.dim_cancer.typing)
  train_X_dim.reduced = z$X_train_reduced
  test_X_dim.reduced = z$X_test_reduced
  cat(sprintf("#reduced dimension: %d\n", n.reduced.dim_cancer.typing))
}
  
model <- fit.ovr.linear.svm(train_X_dim.reduced, train_y, classifier_cancer.typing)
preds <- predict.by.ovr.linear.svm(model, test_X_dim.reduced)
preds$pred.score <- preds$pred.score %>% dplyr::select(all_of(cancer.types)) # reorder columns by the specified order

pred.scores <- data.frame(matrix(nrow = 0, ncol = 3 + length(cancer.types)))
colnames(pred.scores) <- c(cancer.types, "ratio", "pred_label", "sample_id")

# Generate two additional columns: ratio of two largest prediction scores and pred_label
cols_pred.scores <- cancer.types
pred.scores <- preds$pred.score %>%
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
  mutate(pred_label = colnames(.)[apply(.[, cols_pred.scores], 1, which.max)]) %>%  # Identify max column name
  mutate(sample_id = rownames(.))
  
## Add cancer typing prediction columns: prediction scores, ratio and pred_label to "test_samples_for_typing", which are right-adjacent to cancer_type
test_samples_for_typing <- test_samples_for_typing %>%
  left_join(pred.scores, by = "sample_id")

# Output prediction scores and true labels for cancer typing
write.table(test_samples_for_typing, file = out_file, sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat(sprintf("Output: %s\n", out_file))
cat("   pred columns (last 6 columns): liver_cancer, lung_cancer, ovarian_cancer, stomach_cancer, ratio, pred_label\n")
