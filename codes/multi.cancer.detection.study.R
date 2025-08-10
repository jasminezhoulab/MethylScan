library(dplyr)
library(caret)
library(tidyr)
library(stringr)

source("../common_funs.R")

study_id <- "multi.cancer.detection.study"

debug=FALSE
# debug=TRUE

if (debug) {
  input_samples_cv_split_file = ""
  input_data_matrix_file = ""
} else {
  args <- commandArgs(trailingOnly = T)
  
  ## Sample annocation and 5-fold CV split file
  ## Columns: sample_id, true_label (cancer or normal), cancer_type (cancer types or normal), fold_id (Fold1, ..., Fold5)
  input_samples_cv_split_file <- args[1]
  
  ## Data matrix file
  ## Rows are samples and columns are features
  ## Columns: sample_id, feature1, feature2, ...
  input_data_matrix_file <- args[2]
}

classifier = "LinearSVC_l2_c1"
n.reduced.dim = 500

out_dir <- "./output"
dir.create(out_dir, recursive = TRUE)

## Sample annocation and 5-fold CV split file
## Columns: sample_id, true_label (cancer or normal), cancer_type (cancer types or normal), fold_id (Fold1, ..., Fold5)
df.samples.cv <- read.csv(input_samples_cv_split_file)

## Data matrix file
## Rows are samples and columns are features
## Columns: sample_id, feature1, feature2, ...
df.data <- read.csv(input_data_matrix_file, row.names = 1, check.names = FALSE)
df.data[] <- lapply(df.data, function(x) as.numeric(as.character(x)))

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
  
  model <- fit.binary.linear.svm(train_X_dim.reduced, train_y, classifier)
  pred.scores <- predict.by.binary.linear.svm(model, test_X_dim.reduced)
  names(pred.scores) <- test_samples$sample_id
  all.pred.scores <- c(all.pred.scores, pred.scores)
}

df.samples.cv[names(all.pred.scores), "pred_score"] = unname(all.pred.scores)

auc_overall = calculate_auc(label_map[df.samples.cv$true_label],
                                 df.samples.cv$pred_score)
cat(sprintf("AUC_overall: %.4f\n", auc_overall))

write.table(df.samples.cv,
            file = sprintf("%s/multi.cancer.detection.study_pred.csv", out_dir),
            sep = ",",
            row.names = F,
            col.names = T)

cat(sprintf("AUC_overall: %.4f\n", auc_overall),
    file = sprintf("%s/multi.cancer.detection.study_auc.csv", out_dir))

