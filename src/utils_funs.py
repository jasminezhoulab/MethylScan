import sys, re
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

# x is a float numpy array
def sigmoid(x, a=10, b=0):
    return (1 / (1 + np.exp(-a * x + a * b)))


def find_two_largest_values_of_a_list(l):
    ll = sorted(l)  # sort in increasing order
    if len(l) >= 2:
        return ([ll[-1], ll[-2]])
    else:
        return (l)

def get_ratio_of_two_largest_values_of_a_list(l):
    first_largest, second_largest = find_two_largest_values_of_a_list(l)
    if second_largest == 0:
        return np.PINF
    else:
        return (first_largest / float(second_largest))


def apply_dim_reduction(train_df: pd.DataFrame,
                        test_df: pd.DataFrame,
                        dim_arg: str,
                        center_or_standardize: str = "standardize",
                        fold_id: str = ""):
    """
    Apply PCA-based dimensionality reduction.

    Parameters
    ----------
    train_df, test_df : pd.DataFrame
        Rows = samples, columns = features.

    dim_arg : str
        Example: "svd100" or "svd50". The numeric value specifies the
        target dimension (number of principal components).
        If "none", "no", or empty, PCA is skipped.

    fold_id : str
        Optional label for print messages (e.g., "Fold1").

    center_or_standardize : str
        Either "center" or "standardize".
        - "center": subtract mean per feature only (covariance-based PCA)
        - "standardize": subtract mean and divide by std per feature
          (correlation-based PCA, typical z-score)

    Returns
    -------
    train_reduced : pd.DataFrame
        PCA-transformed training data.
    test_reduced : pd.DataFrame
        PCA-transformed test data.
    n_components_used : int
        The actual number of principal components used.

    Notes
    -----
    - PCA always centers internally, but scaling to unit variance is
      performed only if center_or_standardize="standardize".
    - The number of components is capped by both #samples and #features.
    """
    if dim_arg is None:
        dim_arg = "none"
    dim_arg = dim_arg.strip().lower()

    if dim_arg in ("none", "na", "no", ""):
        return train_df, test_df, 0

    # Parse e.g. "svd100" -> 100
    m = re.fullmatch(r"svd(\d+)", dim_arg)
    if not m:
        raise ValueError(f"Invalid dim_reduction arg: {dim_arg}. Use 'svd100', 'svd50', or 'none'.")

    k_req = int(m.group(1))
    n_train, n_feat = train_df.shape
    n_components = min(k_req, n_train, n_feat)

    if n_components < k_req:
        print(f"[{fold_id}] dim_reduction={dim_arg}: "
              f"requested {k_req}, limited to min(#train={n_train}, #features={n_feat}) "
              f"→ using {n_components} components.")

    if n_components <= 0:
        print(f"[{fold_id}] PCA skipped (n_components={n_components}).")
        return train_df, test_df, 0

    # ---- Preprocessing choice ----
    if center_or_standardize.lower() == "standardize":
        scaler = StandardScaler(with_mean=True, with_std=True)
        train_proc = scaler.fit_transform(train_df)
        test_proc = scaler.transform(test_df)
        print(f"[{fold_id}] Standardizing features (mean=0, std=1) before PCA.")
    elif center_or_standardize.lower() == "center":
        # Center manually, since PCA itself centers by default, but
        # we apply the same transform to test data explicitly.
        means = train_df.mean(axis=0)
        train_proc = train_df - means
        test_proc = test_df - means
        print(f"[{fold_id}] Centering features (mean=0) before PCA.")
    elif center_or_standardize.lower() == "none":
        train_proc = train_df
        test_proc = test_df
    else:
        raise ValueError("center_or_standardize must be either 'center', 'standardize', or 'none'.")

    # ---- PCA ----
    pca = PCA(n_components=n_components, svd_solver="auto", random_state=0)
    Xtr = pca.fit_transform(train_proc)
    Xte = pca.transform(test_proc)

    cols = [f"PC{i+1}" for i in range(n_components)]
    train_red = pd.DataFrame(Xtr, index=train_df.index, columns=cols)
    test_red = pd.DataFrame(Xte, index=test_df.index, columns=cols)
    return train_red, test_red, n_components



####################################################################################
####################################################################################
## Multi-class classification for cancer typing
####################################################################################
####################################################################################

def make_multiclass_model(method,
                          nfeature=None,
                          sample_weight=None,
                          num_cpu_used=1,
                          verbose=0):
    items = method.split('_')
    if 'ovr_LinearSVC_' in method:
        # 'ovr_LinearSVC_l2_c1'
        penalty = items[2]
        if penalty == 'l1':
            # dual=True is not implemented for penalty='l1' yet.
            dual = False
        else:
            dual = True
        c = float(items[3][1:])
        if verbose > 0:
            print('  OvR linear SVM (penalty=%s, c=%g)' % (penalty, c))
        model = LinearSVC(penalty=penalty,
                          C=c,
                          dual=dual,
                          multi_class='ovr',
                          class_weight='balanced',
                          max_iter=10000000,
                          random_state=0,
                          verbose=0)
    else:
        sys.stderr.write('Error (make_multiclass_model): options %s is incorrect!\nExit.\n' % method)
        sys.exit(-1)
    return(model)

def classify_multiclass(method,
                        train_x,
                        train_y,
                        test_x=None,
                        num_cpu_used=-1,
                        verbose=0,
                        output_model=False,
                        sample_weight=None,
                        class_order=None    # usually used for cross validation to enforce class order of prediction scores consistent btw different folds
                        ):
    nsample, nfeature = train_x.shape
    model = make_multiclass_model(method,
                                nfeature,
                                sample_weight,
                                num_cpu_used,
                                verbose)

    if sample_weight is None:
        model.fit(train_x, train_y)
    else:
        model.fit(train_x, train_y, sample_weight=sample_weight)

    if test_x is not None:
        if 'SVC' in method:
            # linearSVC do not support predict_proba. decision_function provides the signed distance of that sample to the hyperplane.
            pred = model.decision_function(test_x)
        else:
            pred = model.predict_proba(test_x)
        
        # ---- Ensure consistent class order in score outputs ----
        # Per-class score matrix: align to class_order if provided
        if hasattr(model, "classes_") and pred is not None and pred.ndim == 2:
            model_classes = list(model.classes_)
            if class_order is not None:
                # verify shape matches #classes
                if pred.shape[1] != len(model_classes):
                    raise ValueError(f"Score matrix width ({pred.shape[1]}) "
                                     f"!= #model classes ({len(model_classes)}).")
                try:
                    reindex = [model_classes.index(c) for c in class_order]
                except ValueError as e:
                    raise ValueError(f"class_order contains a class not seen by the model: {e}")
                pred = pred[:, reindex]  # reorder columns to desired order
        
        if output_model:
            return( (pred, model) )
        else:
            return(pred)
    else:
        return(model)
