library(dplyr)
library(rlang)

source("../src/utility_funs.R")

study_id <- "multi.cancer.study_detection"

args <- commandArgs(trailingOnly = TRUE)

## Sample annocation and 5-fold CV split file
## Columns: sample_id, cancer_status (cancer or normal), column_test.or.train (value is a string that indicates test samples/rows)
input_samples_split_file <- args[1]

# Column of the test or train set info. For exmaple, an input string "train_or_test:test" indicates all samples (rows) with value "test" of column name "train_or_test" in "input_samples_split_file" are testing samples; and the rest are training samples
column_test.or.train_info <- args[2] 
parts <- strsplit(column_test.or.train_info, ":")[[1]]
column_test.or.train <- parts[1]
test_value <- parts[2]

## Feature matrix file for multi-cancer detection
## Rows are samples and columns are cancer-specific features
## Columns: sample, feature1, feature2, ...
input_feature_matrix_file_for_cancer.detection <- args[3]

## Output prediction file name
out_file <- args[4]

cat(c(sprintf("Study: %s", study_id),
      sprintf("Sample split file: %s", input_samples_split_file),
      sprintf("Column (test or train indicator): %s (%s)", column_test.or.train, test_value)),
      sprintf("Feature matrix file (row as sample; column as feature): %s", input_feature_matrix_file_for_cancer.detection),
      sep = "\n")

classifier_cancer.detection = "LinearSVC_l2_c1" # Binary Linear SVM with L2 regularization and Cost=1
n.reduced.dim_cancer.detection = 500

## Sample annotations and train/test split file
## Columns: sample_id, column_cancer.status (cancer or normal), column_test.or.train (train or test)
df.samples <- read.csv(input_samples_split_file, header=TRUE)
rownames(df.samples) <- df.samples$sample_id

########################################################################
## Perform multi-cancer detection

## Feature matrix file for multi-cancer detection
## Rows are samples and columns are pancancer features
## Columns: sample_id, feature1, feature2, ...
if (grepl("\\.gz$", input_feature_matrix_file_for_cancer.detection, ignore.case = TRUE)) {
  df.data <- read.csv(gzfile(input_feature_matrix_file_for_cancer.detection), row.names = 1, check.names = FALSE)
} else {
  df.data <- read.csv(input_feature_matrix_file_for_cancer.detection, row.names = 1, check.names = FALSE)
}
df.data[] <- lapply(df.data, function(x) as.numeric(as.character(x)))

cat(c(sprintf("#markers: %d", ncol(df.data)),
      sprintf("Classifier: %s", classifier_cancer.detection)),
      sep = "\n")

label_map <- setNames(c(0, 1), c("normal", "cancer"))

train_samples <- df.samples %>% dplyr::filter(!!sym(column_test.or.train) != test_value)
train_X <- df.data[train_samples$sample_id, ]
train_y <- label_map[train_samples$cancer_status]
cat(sprintf("#train_set: %d (+%d, -%d)\n", nrow(train_X), sum(train_y==1), sum(train_y==0)))

test_samples <- df.samples %>% dplyr::filter(!!sym(column_test.or.train) == test_value)
test_X <- df.data[test_samples$sample_id, ]
test_y <- label_map[test_samples$cancer_status]
cat(sprintf("#test_set: %d (+%d, -%d)\n", nrow(test_X), sum(test_y==1), sum(test_y==0)))

if (n.reduced.dim_cancer.detection >= ncol(train_X)) {
  train_X_dim.reduced = train_X
  test_X_dim.reduced = test_X
  cat(sprintf("DO NOT perform dimension reduction, because #marker=%d <= reduced_dimension=%d\n", ncol(train_X), n.reduced.dim_cancer.detection))
} else {
  z = svd_train_test(train_X, test_X, n.reduced.dim_cancer.detection)
  train_X_dim.reduced = z$X_train_reduced
  test_X_dim.reduced = z$X_test_reduced
  cat(sprintf("#reduced dimension: %d\n", n.reduced.dim_cancer.detection))
}

  
model <- fit.binary.linear.svm(train_X_dim.reduced, train_y, classifier_cancer.detection)
pred.scores <- predict.by.binary.linear.svm(model, test_X_dim.reduced)
names(pred.scores) <- test_samples$sample_id

test_samples[names(pred.scores), "pred.score_multi.cancer"] = unname(pred.scores)

# Output prediction scores for cancer detection
write.table(test_samples, file = out_file, sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat(sprintf("Output: %s\n", out_file))
cat("   pred column (last one column): pred.score_multi.cancer\n")
