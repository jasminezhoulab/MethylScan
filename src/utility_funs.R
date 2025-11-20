library(dplyr)
library(caret)
library(tidyr)
library(stringr)
library(rsample)
library(limma)
library(matrixStats)

# Function to remove rows with all NA values
remove_all_na_rows <- function(df) {
  df[!apply(df, 1, function(row) all(is.na(row))), ]
}

# Sigmoid function
sigmoid <- function(x, a = 3, b = 0) {
  1 / (1 + exp(-a * x + a * b))
}

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

  ## ---- Step to detect sign orientation so that class 1 -> positive by using training data ----
  ## The real reasons why the LiblineaR Linear SVM sometimes gives “flipped” scores — predicting class-0 samples with large positive values and class-1 samples with large negative values — even when the training data are perfectly separable. This is expected behavior from LiblineaR, and not a bug in your code or data. Why scores flip in Linear SVM (LiblineaR):
  ## 1. LiblineaR does not preserve your class label order (0 vs 1). LiblineaR internally sorts class labels and assigns them to its own internal ordering, which may flip the sign of the decision function.
  ## 2. Decision value is not a probability — it is raw distance to hyperplane. The sign depends entirely on how Liblinear sets the direction of w. Even on the same dataset, Liblinear can flip w if: (1) The optimization path converges differently; (2) Class weights (w_i) heavily penalize one class; (3) The cost parameter changes. Since multiplying (w, b) by −1 yields the same classifier, the solver is free to choose either orientation. So flipped signs do NOT mean wrong classification.
  ## 3. Using loss function type = 1 (dual L2-regularized L2-loss SVC) also allows w-sign flipping. Some solvers (e.g., type=1, type=3) have no constraint on sign orientation. So SVM solutions are equivalent under: (w,b) and (−w,−b). Liblinear is allowed to choose either.
  idx_pos_all <- which(y.train == 1)
  if (length(idx_pos_all) > 0) {
    # Use at least 3 training samples (or all if fewer than 3)
    k <- min(3, length(idx_pos_all))
    idx_use <- idx_pos_all[1:k]
    dv_vec <- predict(
      model,
      newx = as.matrix(X.train[idx_use, , drop = FALSE]),
      decisionValues = TRUE
    )$decisionValues[, 1]
    # Average decision values from the 3+ samples
    dv_mean <- mean(dv_vec, na.rm = TRUE)
    model$flip_sign <- !is.na(dv_mean) && dv_mean < 0
  } else {
    model$flip_sign <- FALSE
  }
  model$positive_label <- 1

  return(model)
}

# Function to predict scores using the trained SVM model "LinearSVC_l2_c1"
predict.by.binary.linear.svm <- function(model, X.test) {
  library(LiblineaR)
  # Check if inputs are valid
  # Predict using the model
  prediction <- predict(
    model,
    newx = as.matrix(X.test),
    decisionValues = TRUE
  )
  scores <- prediction$decisionValues[, 1]
  # Extract decision values (scores)
  if (!is.null(model$flip_sign) && isTRUE(model$flip_sign)) {
    scores <- -scores
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

# Function to calculate AUC for a given true_labels and predictions
calculate_auc <- function(true_labels, predictions, quiet = TRUE, direction="<") {
  library(pROC)
  # Check if there are both 1s and 0s in true_labels
  if (length(unique(true_labels)) < 2) {
    return(NA)  # Return NA if only one class is present
  }
  roc_obj <- roc(true_labels, predictions, quiet = quiet, levels = c(0, 1), direction = direction)
  return(as.numeric(auc(roc_obj)))
}

# Function to calculate AUC and its 95% confidence interval
calculate_auc_ci <- function(true_labels, predictions, conf.level=0.95, direction="<") {
  library(pROC)
  # Check if there are both 1s and 0s in true_labels
  if (length(unique(true_labels)) < 2) {
    return(NA)  # Return NA if only one class is present
  }
  roc_obj <- roc(true_labels, predictions, levels = c(0, 1), direction = direction)
  if (!is.na(conf.level)) {
    auc_ci <- ci.auc(roc_obj, conf.level=conf.level)
    return(c(as.numeric(auc(roc_obj)),
             auc_ci[1],
             auc_ci[3]))
  } else (
    return(as.numeric(auc(roc_obj)))
  )
}

#
# Calculate sensitivity and specificity and 95% confidence interval of the sensitivity (optional)
#
# Example usage:
# true_labels <- c(1, 0, 1, 1, 0, 1, 0, 1, 1, 0)
# predictions <- c(0.8, 0.3, 0.6, 0.4, 0.2, 0.8, 0.7, 0.2, 0.3, 0.6)
# max_fp <- 10
# result <- calculate_sens_spec(true_labels, predictions, max_fp, need_original=TRUE, conf.level=0.95)
# print(result)
calculate_sens_spec <- function(true_labels, predictions,
                                max_fp=10, need_original=FALSE,
                                conf.level = NA) {
  library(pROC)
  # Check if there are both 1s and 0s in true_labels
  if (length(unique(true_labels)) < 2) {
    return(NA)  # Return NA if only one class is present
  }
  # Generate ROC object
  roc_obj <- roc(true_labels, predictions, levels = c(0, 1), direction = "<")
  # Get all the ROC coordinates
  all_coords <- coords(roc_obj, "all", ret = c("threshold", "specificity", "sensitivity", "fp", "tp"), transpose = FALSE)
  colnames(all_coords)[4] = "false_positive"
  colnames(all_coords)[5] = "true_positive"
  if (max_fp > max(all_coords$false_positive)) {
    max_fp <- max(all_coords$false_positive)
  }
  # Extract the sensitivities and specificities for given false positive counts
  if (!is.na(conf.level)) {
    sens_spec <- data.frame(
      false_positive = 0:max_fp,
      sensitivity = NA,
      lower_bound_sens = NA,
      upper_bound_sens = NA,
      conf_level = conf.level,
      specificity = NA,
      threshold = NA
    )
  } else {
    sens_spec <- data.frame(
      false_positive = 0:max_fp,
      sensitivity = NA,
      specificity = NA,
      threshold = NA
    )
  }
  rownames(sens_spec) = sens_spec$false_positive
  # Total number of positive cases
  n_positive <- sum(true_labels == 1)
  # Loop through false positive numbers and find the closest match
  for (i in 0:max_fp) {
    # Find the closest row where false positives are approximately equal to `i`
    closest_row <- all_coords[which.min(abs(all_coords$false_positive - i)), ]
    sens_spec[as.character(i), "false_positive"] <- closest_row$false_positive
    sens_spec[as.character(i), "sensitivity"] <- closest_row$sensitivity
    sens_spec[as.character(i), "specificity"] <- closest_row$specificity
    sens_spec[as.character(i), "threshold"] <- closest_row$threshold
    # If conf.level is provided, calculate CI for sensitivity using `true_positive`
    if (!is.na(conf.level)) {
      true_positive <- closest_row$true_positive # True positives from all_coords
      ci <- clopper_pearson_ci(n_positive, true_positive, conf.level)
      sens_spec[as.character(i), "lower_bound_sens"] <- ci[1]
      sens_spec[as.character(i), "upper_bound_sens"] <- ci[2]
    }
  }
  if (need_original) {
    return(list(sens_spec=sens_spec, all=all_coords))
  } else {
    return(sens_spec)
  }
}

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
  cancer_type_col <- sym(cancer_type_colum.name)
  pred_col <- sym(pred_column.name)
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

################################################################# 
############# Marker Discovery

## Select cancer markers by the methylation matrix of the paired tissue samples (tumors and their adjacent normal tissues)
## tumor_tissues_list: list of tumor samples whose sample_id contains a substring "_T", which could be replaced by "_N" that represents sample_id of the adjacent normal tissue sample.
## tissue_matrix: a matrix (or data frame) with rows as markers and columns as tissue samples. It has both rownames and colnames.
## diff_cutoff: the threshold that is required for the methylation level difference between the tumor and its adjacent normal tissue must be > diff_cutoff.
## topK: Choose the top "topK" markers whose occurrence or frequency (number of tissue pairs whose methylation level difference >= diff_cutoff) is one of the largest "topK".
select_markers_by_paired.tissues_and_occurrence <- function(tumor_tissues_list,
                                                            tissue_matrix,
                                                            diff_cutoff,
                                                            topK) {
  candidate_markers=rownames(tissue_matrix)
  candidate_markers = candidate_markers[!grepl("chrX", candidate_markers)]
  candidate_markers = candidate_markers[!grepl("chrY", candidate_markers)]
  candidate_markers = candidate_markers[!grepl("chrM", candidate_markers)]
  tissue_matrix <- tissue_matrix[candidate_markers, ]
  diff_matrix <- data.frame(tissue_matrix[, tumor_tissues_list])
  for (tumor_sample in tumor_tissues_list){
    normal_sample <- gsub("_T", "_N", tumor_sample)
    diff_matrix[,tumor_sample] <- tissue_matrix[,tumor_sample] - tissue_matrix[,normal_sample]
  }
  diff_matrix$num_over_dif_thresh <- apply(diff_matrix[, tumor_tissues_list], 1, function(x) {sum(x > diff_cutoff)})
  topN_th = sort(diff_matrix$num_over_dif_thresh, decreasing = TRUE)[topK]
  if (topN_th < 1) {
    topN_th = 1
  } 
  diff_matrix <- subset(diff_matrix, diff_matrix$num_over_dif_thresh >= topN_th)
  markers <- rownames(diff_matrix)
  return(markers)
}

## Rows are features, columns are samples
# Remove those rows whose number of NA exceeds half of all columns
transform_data_rows <- function(df) {
  sum.na <- function(x){ return(sum(is.na(x))) }
  # This line counts NAs per row (correct)
  feature.na.cnt <- apply(df, 1, sum.na)
  # This line removes rows where ALL values are NA
  df <- df[feature.na.cnt != ncol(df),]
  # This line recalculates NA proportion per row (after dropping all-NA rows)
  data.na.cnt <- apply(df, 1, sum.na) / ncol(df)
  # This line removes rows with >= 50% NAs (this is the core filter)
  return(df[data.na.cnt < 0.5,])
}


## Rows are samples, columns are features
# Remove those columns whose number of NA exceeds half of all rows
transform_data_cols <- function(df) {
  sum.na <- function(x) { sum(is.na(x)) }
  # Step 1: Count NA per column
  feature.na.cnt <- apply(df, 2, sum.na)  # apply over columns
  # Step 2: Remove columns where all values are NA
  df <- df[, feature.na.cnt != nrow(df)]
  # Step 3: Recalculate % NA per column (after dropping all-NA cols)
  data.na.cnt <- apply(df, 2, sum.na) / nrow(df)
  # Step 4: Keep columns with less than 50% missing
  df <- df[, data.na.cnt < 0.5]
  return(df)
}

make.contrasts.btw.two.classes <- function(classes) {
  contr = c()
  for (i in classes) {
    for (j in classes) {
      if (i != j) {
        contr = c(contr, paste0(i, '-', j))
      }
    }
  }
  return(contr)
}

make.contrasts.btw.multi.classes <- function(classes, num.positive.classes) {
  contr = c()
  if (num.positive.classes==1) {
    # only 1 class is plus, all other classes are minus
    for (c in classes) {
      # for example, "adipose-b_cell-colon-esophagus-heart-kidney-liver-lung-monocyte-neutrophils-pancreas-small_intestine-spleen-stomach-t_cell"
      formula_ = paste0(c,paste0('-',classes[!(classes==c)],collapse=''),sep='')
      contr = c(contr, formula_)
    }
  } else if (num.positive.classes>=length(classes)) {
    stop('num.positive.classes must be < length(classes)')
  } else {
    # columns are combinations of classes
    comb.classes = combn(classes, num.positive.classes)
    for (j in 1:ncol(comb.classes)) {
      # the combination of 'num.positive.classes' classes is plus, all other classes are minus
      plus.classes = paste0(comb.classes[,j],collapse='+')
      minus.classes = paste0('-',setdiff(classes,comb.classes[,j]),collapse='')
      contr = c(contr, paste0(plus.classes, minus.classes))
    }
  }
  return(contr)
}

select_top_general <- function(df, cutoff=0.7, num=50) {
  filter_df <- df[abs(df$class.diff) > cutoff,]
  filter_df <- filter_df[!grepl("^chrX|^chrY|^chrM", filter_df$feature.list), ]
  df_sorted <- filter_df[order(abs(filter_df$logFC), decreasing = TRUE),]
  return(head(df_sorted, num))
}

##
## Use Limma method to select one-vs-rest features for multi-class data
##   annot: a data frame of training samples, with three major columns
##      (1) sample: training samples names
##      (2) type: class names
##   data: a data frame with rows as features and columns as samples, which is "data.beta.imp" in the original code
##
select_ovr_features <- function(annot, data,
                                cutoff=0.7,
                                topk.features=50,
                                output_vector=TRUE,
                                quiet=TRUE) {
  classes_list <- sort(unique(annot$type))
  n_class <- length(classes_list)
  
  design <- model.matrix(~0+type, data=annot)
  # fit the linear model
  vfit <- lmFit(data, design, na.action=na.omit)
  
  # remove NA for the linear model
  coef <- vfit$coefficients
  sum.na <- function(x){return(sum(is.na(x)))}
  coef.na.cnt <- apply(coef, 1, sum.na)
  valid.marker.id <- names(which(coef.na.cnt == 0))
  data <- data[valid.marker.id, ]
  
  # fit the linear model again
  vfit <- lmFit(data, design, na.action=na.omit)
  
  # create a contrast matrix for specific comparisons
  contrasts <- c()
  for (num.positive.class in c(1)) {
    contrasts = c(contrasts, make.contrasts.btw.multi.classes(classes_list, num.positive.class))
  }
  contr.matrix <- makeContrasts(contrasts=contrasts, levels=classes_list)
  # fit the contrasts
  vfit <- contrasts.fit(vfit, contr.matrix)
  tfit <- treat(vfit, lfc = 1, trend = F, robust = F)
  results.bh <- decideTests(tfit)
  
  MARKER.LIST = list()
  for(index_class in 1:length(contrasts)){
    topK = dim(data)[1]
    t = classes_list[index_class]
    
    if (n_class > 2) {
      treat.table <- topTreat(tfit,num=topK,coef=index_class,genelist=rownames(data),sort.by="p", resort.by="logFC")
    } else if (n_class == 2) {
      treat.table <- topTreat(tfit,num=Inf,coef=index_class,genelist=rownames(data),sort.by="p", resort.by="logFC")
      treat.table <- treat.table[treat.table$logFC > 0, ]
      topK <- nrow(treat.table)
    }
    
    a = rownames(treat.table)
    logFC <- treat.table$logFC
    class.mean <- apply(as.matrix(data[a,t==annot[,'type']]), 1, mean, na.rm = T)
    others.mean <- apply(as.matrix(data[a,t!=annot[,'type']]), 1, mean, na.rm = T)
    class.na <- apply(is.na(as.matrix(data[a,t==annot[,'type']])), 1, sum)
    others.na <- apply(is.na(as.matrix(data[a,t!=annot[,'type']])), 1, sum)
    feature.list <- a
    class.diff <- abs(class.mean - others.mean)
    t.stats <- treat.table$t
    adj.p.values <- treat.table$adj.P.Val
    dt <- data.frame(classes_list = rep(t, each = topK),
                     feature.list = feature.list, class.diff = class.diff,
                     class.mean = class.mean, others.mean = others.mean,
                     class.na = class.na, others.na = others.na,
                     t.stats = t.stats, adj.p.values = adj.p.values,
                     logFC = logFC, stringsAsFactors = F)
    MARKER.LIST[[t]] <- dt
  }

  # Initialize an empty list to store the marker vectors
  marker_list <- list()
  
  # Iterate over each name in the MARKER.LIST
  for (class_name in names(MARKER.LIST)) {
    marker_list[[class_name]] <- rownames(select_top_general(MARKER.LIST[[class_name]], cutoff, topk.features))
  }
  
  # Create the summary data frame with lengths
  z <- data.frame(
    t(sapply(marker_list, length))
  )
  if (!quiet) {print(z)}
  
  # Combine all markers into a single unique vector
  union_ovr <- unique(unlist(marker_list))
  
  if (output_vector) {
    return(union_ovr)
  } else {
    return(list(union=union_ovr, stat=z, details=marker_list))
  }
}
