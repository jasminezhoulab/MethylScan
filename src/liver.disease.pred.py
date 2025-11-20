import os, argparse, csv
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from utils_funs import sigmoid, apply_dim_reduction, classify_multiclass, get_ratio_of_two_largest_values_of_a_list

def main(debug=False):
    if debug:
        # Set the working directory to the directory of the current Python file
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        
        # Provide inline arguments for debugging
        cv_splits_file = "../samples/CV.5folds/CV.5folds.by.etiology.age.seed345_Repeat01.csv"
        sample_column = "plasma_sample"
        class_column = "liver_disease"
        column_test_or_train_info = "train_or_test:test"  # "<column>:<test_value>"

        classifier = "ovr_LinearSVC_l2_c1"
    
        feature_file = "../markers/CV.5folds/ovo.markers.plasma.limma.topk50.cut0.9_CV.5folds.by.etiology.age.seed345_Repeat01.[fold_id].csv"
        
        dim_reduction = "svd100"  # or "svd100" / "svd50" / "none" when testing
        center_or_standardize_before_dim_reduction = "standardize"
            
        input_data_file = f'../../Wenyuan.replicate.Ran.work/data/output/liver.log1k.data.csv.gz'
    
        out_pred_file = './output/pred.csv'
    else:
        parser = argparse.ArgumentParser(description='Multi-class classification')
        parser.add_argument('--classifier', type=str, required=True, help='Classification Method')
        parser.add_argument('--feature_file', type=str, required=True, help='Input feature file (one feature per line)')
        parser.add_argument('--dim_reduction', type=str, default='none',
                    help="Dimension reduction option: 'svd100', 'svd50', or 'none' (default).")
        parser.add_argument('--center_or_standardize_before_dim_reduction', type=str, default='standardize',
                            choices=['center', 'standardize', 'none'], help="Preprocess features before PCA.")
        parser.add_argument('--sample_split_file', type=str, required=True, help='Input sample split (train and test set) csv file with clinical info')
        parser.add_argument('--sample_column', type=str, required=True, help='Sample column name in the input clinical and cross-validation csv file')
        parser.add_argument('--class_column', type=str, required=True, help='Sample class column name in the input clinical and cross-validation csv file')
        parser.add_argument('--column_test_or_train_info', type=str, required=True, help='Format "<column>:<test_value>" (e.g., "train_or_test:test"). Column whose value "test" marks test samples; others are training')
        parser.add_argument('--input_data_file', type=str, required=True, help='Input data csv file (rows are samples, columns are features)')
        
        parser.add_argument('--out_pred_file', type=str, required=True, help='Output predictions CSV file')

        args = parser.parse_args()
        classifier = args.classifier
        feature_file = args.feature_file
        dim_reduction = args.dim_reduction
        center_or_standardize_before_dim_reduction = args.center_or_standardize_before_dim_reduction
        sample_split_file = args.sample_split_file
        sample_column = args.sample_column
        class_column = args.class_column
        column_test_or_train_info = args.column_test_or_train_info
        input_data_file = args.input_data_file
        out_pred_file = args.out_pred_file

    print(f"{'-'*40}")
    print(f"{'Predict liver disease: Parameters':^40}")
    print(f"{'-'*40}")
    print(f"{'classifier':<28}: {classifier}")
    print(f"{'feature_file':<28}: {feature_file}")    
    print(f"{'dim_reduction':<28}: {dim_reduction}")
    if dim_reduction != "none":
        print(f"{'center_or_standardize_before_dim_reduction':<28}: {center_or_standardize_before_dim_reduction}")
    print(f"{'sample_split_file':<28}: {sample_split_file}")
    print(f"{'  sample_column':<28}: {sample_column}")
    print(f"{'  class_column':<28}: {class_column}")
    print(f"{'  column_test_or_train_info':<28}: {column_test_or_train_info}")
    print(f"{'input_data_file':<28}: {input_data_file}")
    print(f"{'out_pred_file':<28}: {out_pred_file}")
    print(f"{'-'*40}")
    
    # --- Parse "<column>:<test_value>" ---
    parts = column_test_or_train_info.split(":")
    if len(parts) != 2 or not parts[0] or not parts[1]:
        raise ValueError(
            f'--column_test_or_train_info must be in the form "<column>:<test_value>", '
            f'e.g. "train_or_test:test"; got "{column_test_or_train_info}"'
        )
    column_test_or_train = parts[0]
    test_value = parts[1]

    # Load sample split and clinical data
    annot = pd.read_csv(sample_split_file, dtype=str)
    if column_test_or_train not in annot.columns:
        raise ValueError(f'Column "{column_test_or_train}" not found in {sample_split_file}')
    if sample_column not in annot.columns:
        raise ValueError(f'Column "{sample_column}" not found in {sample_split_file}')
    if class_column not in annot.columns:
        raise ValueError(f'Column "{class_column}" not found in {sample_split_file}')
    
    annot.index = annot[sample_column]
    annot["type"] = annot[class_column]
    
    # Determine train/test sets (case-insensitive match on test_value)
    is_test = (annot[column_test_or_train].astype(str).str.lower()
               == str(test_value).lower())
    test_samples = annot.index[is_test].tolist()
    train_samples = annot.index[~is_test].tolist()
    
    if len(test_samples) == 0:
        raise ValueError(f'No test samples found (no rows where {column_test_or_train} == "{test_value}").')
    if len(train_samples) == 0:
        raise ValueError('No training samples found (all rows matched test).')
    
    # --- Stable class ordering across the whole split file ---
    classes_fixed = sorted(annot["type"].unique())
    le_global = LabelEncoder().fit(classes_fixed)
    class_names = list(le_global.classes_)

    # --- Load the feature matrix (samples x features) ---
    df_methyl_all_data = pd.read_csv(input_data_file, index_col=0)
    df_methyl_all_data = df_methyl_all_data.astype(float, errors='raise')

    # --- Read selected features ---
    with open(feature_file, 'r') as f:
        features = [line.strip() for line in f if line.strip()]
    missing_feats = [ft for ft in features if ft not in df_methyl_all_data.columns]
    if missing_feats:
        raise ValueError(f"{len(missing_feats)} selected features not in data, e.g. {missing_feats[:5]}")

    # --- Subset & align ---
    train_data = df_methyl_all_data.loc[train_samples, features]
    test_data = df_methyl_all_data.loc[test_samples, features]
    
    # --- Optional PCA ---
    train_data, test_data, n_used = apply_dim_reduction(
        train_data, test_data,
        dim_arg=dim_reduction,
        center_or_standardize=center_or_standardize_before_dim_reduction,
        fold_id="single"
    )
    if n_used > 0:
        print(f"[single] PCA applied: n_components={n_used}, "
              f"train_shape={train_data.shape}, test_shape={test_data.shape}")
        
    # --- Labels ---
    train_labels = le_global.transform(annot.loc[train_data.index, "type"])
    test_labels  = le_global.transform(annot.loc[test_data.index, "type"])
        
    # --- Train & predict ---
    predictions = classify_multiclass(
        classifier,
        train_data,
        train_labels,
        test_data,
        class_order=le_global.transform(class_names).tolist()  # fixed column order
    )
        
    if "SVC" in classifier:
        pred_scores_sigmoid = sigmoid(predictions, 3)  # transform scores by sigmoid function, so that all scores have distribution: 0 -> 0.5, +infinity -> 1, -infinity -> 0
        # pred_scores_sigmoid = pred_scores_sigmoid / sum(pred_scores_sigmoid)  # transform scores to probabilities, such that sum of all scores equals to 1.
        
        # row-wise normalize so each sample’s scores sum to 1
        row_sums = pred_scores_sigmoid.sum(axis=1, keepdims=True)
        # avoid divide-by-zero if any row is all zeros
        row_sums[row_sums == 0] = 1.0
        predictions = pred_scores_sigmoid / row_sums 
        
    # --- Build predictions dataframe (no performance metrics) ---
    df_pred = pd.DataFrame(predictions, columns=class_names, index=test_data.index).reset_index()
    df_pred.rename(columns={"index": "sample"}, inplace=True)
    df_pred.insert(1, 'true_label', [class_names[l] for l in test_labels])
    df_pred.insert(2, 'pred_label', df_pred[class_names].idxmax(axis=1))
    # ratio of top-2 scores per row
    df_pred.insert(3, 'ratio',
                   df_pred[class_names].apply(get_ratio_of_two_largest_values_of_a_list, axis=1))

    # --- Write predictions CSV ---
    df_pred.to_csv(out_pred_file, index=False)  
        

if __name__ == "__main__":
    main(debug=False)  # Set to False to use command-line arguments
