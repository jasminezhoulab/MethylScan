library(dplyr)
library(caret)
library(tidyr)
library(stringr)
library(rsample)

################################################################# 
############# Dimension reduction

filter.out.near.zero.variance.columns <- function(X) {
  library(caret)
  # Remove zero-variance or near-zero variance features
  nzv <- caret::nearZeroVar(X,
                            freqCut = 899/1,   # Flag if the most common value is 99x more frequent than the second
                            uniqueCut = 1,     # Flag if unique values < 1% of total samples
                            names = TRUE)
  features.used <- setdiff(colnames(X), nzv)
  return(features.used)
}

# SVD-based Dimensionality Reduction (Training + Testing)
svd_train_test <- function(X_train, X_test, k, verbose=F) {
  library(caret)
  # Remove zero-variance or near-zero variance features
  nzv <- caret::nearZeroVar(X_train,
                            freqCut = 899/1,   # Flag if the most common value is 99x more frequent than the second
                            uniqueCut = 1,     # Flag if unique values < 1% of total samples
                            names = TRUE)
  features.used <- setdiff(colnames(X_train), nzv)
  if (verbose) {
    cat(sprintf("#dim = %d = %d-%d\n", length(features.used), nrow(X_train), length(ncv)))
  }
  X_train <- X_train[, features.used]
  X_test <- X_test[, features.used]
  
  # Step 1: Standardize the training data
  train_means <- colMeans(X_train)  # Compute column means
  train_sds <- apply(X_train, 2, sd)  # Compute column standard deviations
  
  X_train_scaled <- scale(X_train, center = train_means, scale = train_sds)  # Standardize train data
  
  # Step 2: Compute SVD on training data
  svd_result <- svd(X_train_scaled)
  
  # Determine the maximum number of components (ncomp) based on the dataset and predefined ncomp
  k <- min(min(nrow(X_train), ncol(X_train)), k)
  
  # Step 3: Select top k singular vectors (principal directions)
  V_k <- svd_result$v[, 1:k]  # Right singular vectors corresponding to top k singular values
  
  # Step 4: Project training data onto the top k dimensions
  X_train_reduced <- X_train_scaled %*% V_k
  
  # Step 5: Apply the same transformation to test data
  X_test_scaled <- scale(X_test, center = train_means, scale = train_sds)  # Standardize test data
  X_test_reduced <- X_test_scaled %*% V_k  # Project test data onto k dimensions
  
  return(list(X_train_reduced = X_train_reduced, 
              X_test_reduced = X_test_reduced,
              k = k,
              V_k = V_k,
              features.used = features.used,
              train_means = train_means, 
              train_sds = train_sds))
}

################################################################# 
############# Binary classification by Linear SVM

## Binary classification by Linear SVM
fit.binary.linear.svm <- function(X.train, y.train, classifier = "LinearSVC_l2_c1") {
  library(LiblineaR)
  # Check if inputs are valid
  if (!is.numeric(y.train) || length(unique(y.train)) != 2) {
    stop("y.train must be a numeric vector with two unique values (0 and 1).")
  }
  if (nrow(X.train) != length(y.train)) {
    stop("Number of rows in X.train must match the length of y.train.")
  }
  # Parse the classifier string
  classifier_parts <- strsplit(classifier, "_")[[1]]
  if (length(classifier_parts) != 3 || classifier_parts[1] != "LinearSVC") {
    stop("Invalid classifier format. Expected format: 'LinearSVC_l2_cX', where X is the cost value.")
  }
  # Determine the type based on loss function
  loss_function <- classifier_parts[2]
  if (loss_function == "l2") {
    ## L2-regularized L2-loss support vector classification (dual)
    type <- 1
  } else if (loss == "l1") {
    ## L1-regularized L2-loss support vector classification
    type = 5
  } else {
    stop("Error: Unsupported loss function. Classifier (", classifier, ") for LinearSVC must use loss function l1 or l2!\n", sep="")
  }
  # Extract cost value
  cost <- as.numeric(sub("c", "", classifier_parts[3]))
  if (is.na(cost)) {
    stop("Invalid cost value in classifier string.")
  }
  epsilon <- 1e-5 ## default value
  bias <- 1 ## default value for intercept
  
  ## class_weight must be "balanced". Implemented using strategy of https://scikit-learn.org/stable/modules/generated/sklearn.svm.LinearSVC.html
  ## n_samples / (n_classes * np.bincount(y))
  wi <- nrow(X.train) / (2 * count.frequncy(y.train))
  
  # Train the model with the specified parameters
  model <- LiblineaR(
    data = as.matrix(X.train),
    target = as.numeric(y.train),
    type = type,
    cost = cost,
    epsilon = epsilon,
    bias = bias,
    wi = wi,
    verbose = FALSE
  )
  return(model)
}

# Function to predict scores using the trained SVM model "LinearSVC_l2_c1"
# A bug in this LibLineaR, sometimes the trained model can lead to reversed predictions,
# so we need to use "reverse.decision.values" to manually inverse them.
predict.by.binary.linear.svm <- function(model, X.test, reverse.decision.values = F) {
  library(LiblineaR)
  # Check if inputs are valid
  # Predict using the model
  prediction <- predict(
    model,
    newx = as.matrix(X.test),
    decisionValues = TRUE
  )
  # Extract decision values (scores)
  if (reverse.decision.values) {
    scores <- 0 - prediction$decisionValues[, 1]
  } else {
    scores <- prediction$decisionValues[, 1]
  }
  
  # print(prediction$decisionValues)
  return(scores)
}

################################################################# 
############# Multi-class classification by Linear SVM

## count.frequncy(c("hello","hello","hello","hello","hello","world", "world"))
## return a vector:
# hello world
# 5     2
count.frequncy <- function(vec) {
  counts <- table(vec)
  ret <- setNames(as.vector(counts), as.character(names(counts)))
  return(ret)
}

## Multi-class classification by One-vs-Rest Linear SVM: "ovr_LinearSVC_l2_c1"
fit.ovr.linear.svm <- function(X.train, y.train, classifier="ovr_LinearSVC_l2_c1") {
  library(LiblineaR)
  items <- strsplit(classifier, "_")[[1]]
  if (items[2] != "LinearSVC") {
    stop("Error: classifier (", classifier, ") is not LinearSVC!\n", sep="")
  }
  multi_class <- items[1]
  if (multi_class != "ovr") {
    stop("Error: classifier (", classifier, ") for LinearSVC supports only One-vs-Rest (ovr)!\n", sep="")
  }
  loss <- items[3]
  if (loss == "l2") {
    ## L2-regularized L2-loss support vector classification (dual)
    type = 1
  } else if (loss == "l1") {
    ## L1-regularized L2-loss support vector classification
    type = 5
  } else {
    stop("Error: classifier (", classifier, ") for LinearSVC must use loss function l1 or l2!\n", sep="")
  }
  cost <- as.numeric(substr(items[4], 2, nchar(items[4])))
  epsilon <- 1e-4 ## default value
  bias <- 1 ## default value for intercept
  
  n.train.sample <- nrow(X.train)
  n.feature <- ncol(X.train)
  n.class <- length(unique(y.train))
  
  ## class_weight must be "balanced". Implemented using strategy of https://scikit-learn.org/stable/modules/generated/sklearn.svm.LinearSVC.html
  ## n_samples / (n_classes * np.bincount(y))
  wi <- n.train.sample / (n.class * count.frequncy(y.train))
  model <- LiblineaR(X.train,
                     y.train,
                     type = type,
                     cost = cost,
                     epsilon = epsilon,
                     bias = bias,
                     wi = wi,
                     verbose = FALSE)
  return(model)
}

## Make predictions when using "ovr_LinearSVC_l2_c1"
predict.by.ovr.linear.svm <- function(model, X.test) {
  library(LiblineaR)
  ret <- predict(model,
                 X.test,
                 decisionValues = TRUE)
  ret$decisionValues <- data.frame(ret$decisionValues)
  rownames(ret$decisionValues) <- rownames(X.test)
  names(ret$predictions) <- rownames(X.test)
  return( list(pred.score=ret$decisionValues,
               pred.label=ret$predictions))
}

current.time <- function() {
  return(format(Sys.time(), "%Y.%m.%d_%H.%M.%S"))
}

################################################################# 
############# Performance of cancer detection (binary classification)

#
# preds: a data frame with true label column (cancer_type_colum.name) + prediction score column (pred_column.name)
# cancer_type_colum.name: column name (string) that contains each cancer type's name and "normal" class
# pred_column.name: column name (string) of prediction scores
# Return a data frame (ret): 2 row and columns. Column 1 is specificity, Column 2 is threshold, and the following each column is a cancer type; Row 1 and 2 values are AUC value (a string of AUC, a string of AUC (lower_bound, upper_bound) if conf.level is given) and sample size of this cancer type.
#
calc.auc.sens.spec_cancer.types_and_generate_report_for_cancer.detection <- function(
    preds,
    cancer_type_colum.name,
    pred_column.name,
    conf.level=0.95,
    max_fp_for_specificity = 10, direction="<") {
  if (is.null(preds) || (length(preds) == 1 && is.na(preds))) {
    return(NA)
  }
  library(rlang)
  library(dplyr)
  cancer_type_col <- ensym_wenyuan(cancer_type_colum.name)
  pred_col <- ensym_wenyuan(pred_column.name)
  preds_df <- preds %>%
    mutate(true_label_binary = ifelse(grepl("cancer|leukemia|melanoma|lymphoma", !!cancer_type_col), 1, 0))
  # Calculate AUC for each cancer_type
  cancer_types_list <- setdiff(unique(preds_df[[cancer_type_col]]), c("normal"))
  ret <- setNames(data.frame(matrix(ncol = 3+ length(cancer_types_list), nrow = 2+max_fp_for_specificity+1)),
                  c("threshold", "specificity", "all_cancers", cancer_types_list))
  rownames(ret) = c("auc", "sample_size", paste0("fp", 0:max_fp_for_specificity))
  for (c in setdiff(colnames(ret), c("threshold", "specificity"))) {
    cat(c, "\n")
    if (c == "all_cancers") {
      z <- preds_df
    } else {
      z <- rbind(preds_df %>% filter((!!cancer_type_col==UQ(c)) | (!!cancer_type_col=="normal")))
    }
    auc_detail = calculate_auc_ci(z$true_label_binary, z[[pred_col]], conf.level, direction=direction)
    sens_spec_detail = calculate_sens_spec(z$true_label_binary, z[[pred_col]], max_fp=max_fp_for_specificity, conf.level=conf.level)
    if (is.null(sens_spec_detail) || (length(sens_spec_detail) == 1 && is.na(sens_spec_detail))) {
      return(NA)
    }
    if (!is.na(conf.level)) {
      auc = auc_detail[1]
      auc_lower_bound = auc_detail[2]
      auc_upper_bound = auc_detail[3]
      ret["auc", c] = sprintf("%.4g (%.4g - %.4g)", auc, auc_lower_bound, auc_upper_bound)
      for (k in 1:nrow(sens_spec_detail)) {
        fp = sens_spec_detail[k, "false_positive"]
        sensitivity = sens_spec_detail[k, "sensitivity"]
        lower_bound_sens = sens_spec_detail[k, "lower_bound_sens"]
        upper_bound_sens = sens_spec_detail[k, "upper_bound_sens"]
        specificity = sens_spec_detail[k, "specificity"]
        threshold = sens_spec_detail[k, "threshold"]
        ret[paste0("fp",fp), c] = sprintf("%.4g (%.4g - %.4g)", sensitivity, lower_bound_sens, upper_bound_sens)
        ret[paste0("fp",fp), "specificity"] = sprintf("%.4g", specificity) # Specificity for each cancer stage is the same, but different btw different fp
        ret[paste0("fp",fp), "threshold"] = sprintf("%.4g", threshold) # Threshold for each cancer stage is the same, but different btw different fp
      }
    } else {
      auc = auc_detail
      ret["auc", c] = sprintf("%.4g", auc)
      for (k in 1:nrow(sens_spec_detail)) {
        fp = sens_spec_detail[k, "false_positive"]
        sensitivity = sens_spec_detail[k, "sensitivity"]
        specificity = sens_spec_detail[k, "specificity"]
        threshold = sens_spec_detail[k, "threshold"]
        ret[paste0("fp",fp), c] = sprintf("%.4g", sensitivity)
        ret[paste0("fp",fp), "specificity"] = sprintf("%.4g", specificity) # Specificity for each cancer stage is the same, but different btw different fp
        ret[paste0("fp",fp), "threshold"] = sprintf("%.4g", threshold) # Threshold for each cancer stage is the same, but different btw different fp
      }
    }
    ret["sample_size", c] = sum(z$true_label_binary==1, na.rm = T)
  }
  ret <- remove_all_na_rows(ret)
  return(ret)
}
