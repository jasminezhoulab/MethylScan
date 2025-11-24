library(dplyr)
library(rlang)
library(tibble)
library(caret)
library(tidyr)
library(stringr)

source("../src/utility_funs.R")

study_id <- "race.study"

args <- commandArgs(trailingOnly = T)

## Sample annocation and training/testing sample split file
## Columns: sample_id, race, column_test.or.train (value is a string that indicates test samples/rows)
input_samples_split_file <- args[1]

# Column of the test or train set info. For exmaple, an input string "train_or_test:test" indicates all samples (rows) with value "test" of column name "train_or_test" in "input_samples_split_file" are testing samples; and the rest are training samples
column_test.or.train_info <- args[2] 
parts <- strsplit(column_test.or.train_info, ":")[[1]]
column_test.or.train <- parts[1]
test_value <- parts[2]

## Feature matrix file for race prediction study
## Rows are samples and columns are race features
## Columns: sample_id, feature1, feature2, ...
input_feature_matrix_file <- args[3]

## Race markers file. Each line is a marker id.
input_marker_file <- args[4]

## Output prediction file name
out_pred_file <- args[5]

classifier = "ovr_LinearSVC_l2_c1"

## Sample annotation and train/test split file
## Columns: sample_id, race, column_test.or.train (train or test)
df.samples <- read.csv(input_samples_split_file)
rownames(df.samples) <- df.samples$sample_id

races_list <- sort(unique(df.samples$race))

########################################################################
## Perform race prediction

markers <- readLines(input_marker_file)

## Feature matrix file for race prediction study
## Rows are samples and columns are race features
## Columns: sample_id, feature1, feature2, ...
df.data <- read.csv(input_feature_matrix_file, row.names = 1, check.names = FALSE)
df.data[] <- lapply(df.data, function(x) as.numeric(as.character(x)))
df.data <- df.data[, markers]

cat(c(sprintf("#markers: %d", ncol(df.data)),
      sprintf("Classifier: %s", classifier)),
      sep = "\n")

# all.pred.scores <- data.frame(matrix(nrow = 0, ncol = 3 + length(races_list)))
# colnames(all.pred.scores) <- c(races_list, "ratio", "pred.label", "sample_id")

train_samples <- df.samples %>% dplyr::filter(!!sym(column_test.or.train) != test_value)
train_X <- df.data[train_samples$sample_id, ]
train_y <- train_samples$race
counts <- sapply(races_list, \(x) sum(train_samples$race == x))
cat(sprintf("#train_set: %d (%s)\n", nrow(train_X), paste(sprintf("%s: %d", names(counts), counts), collapse = ", ")))

test_samples <- df.samples %>% dplyr::filter(!!sym(column_test.or.train) == test_value)
test_X <- df.data[test_samples$sample_id, ]
test_y <- test_samples$race
counts <- sapply(races_list, \(x) sum(test_samples$race == x))
cat(sprintf("#test_set: %d (%s)\n", nrow(test_X), paste(sprintf("%s: %d", names(counts), counts), collapse = ", ")))

model <- fit.ovr.linear.svm(train_X, train_y, classifier)
preds <- predict.by.ovr.linear.svm(model, test_X)
preds$pred.score <- preds$pred.score %>% dplyr::select(all_of(races_list))

pred.scores <- data.frame(matrix(nrow = 0, ncol = 3 + length(races_list)))
colnames(pred.scores) <- c(races_list, "ratio", "pred_label", "sample_id")

# Generate two additional columns: ratio of two largest prediction scores and pred_label
cols_pred.scores <- races_list
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
  mutate(pred.label = colnames(.)[apply(.[, cols_pred.scores], 1, which.max)]) %>%  # Identify max column name
  mutate(sample_id = rownames(.))

## Add race prediction columns: prediction scores, ratio and pred_label to "test_samples", which are right-adjacent to race
test_samples <- test_samples %>%
  left_join(pred.scores, by = "sample_id")

# Output prediction scores, predicted labels, and true labels
write.table(test_samples, file = out_pred_file, sep = ",", row.names = F, col.names = T, quote = FALSE, )

cat(sprintf("Output: %s\n", out_pred_file))
cat("   pred columns (last multiple columns): race_1, race_2, ..., ratio, pred_label\n")

