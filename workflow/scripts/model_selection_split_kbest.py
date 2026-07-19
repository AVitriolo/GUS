import numpy as np
import pandas as pd
import sklearn as sk
import skopt as sko
import shap as sh
import joblib as jbl
import argparse
import time
import sys

from sklearn.feature_selection import SelectKBest, f_regression

from workflow.scripts.helpers.BUILD_MODEL import *
from workflow.scripts.helpers.GET_HYPERPARAM_SPACE_KBEST import *
from workflow.scripts.helpers.OPTIMIZE_HYPERPARAMS import *
from workflow.scripts.helpers.GET_EXPLAINER import *
from workflow.scripts.helpers.EVALUATE_PERFORMANCE_nocv import *

# Split happens upstream (split_samples.py); do NOT re-split here.

sys.stdout = sys.stderr

parser = argparse.ArgumentParser()
parser.add_argument('--input_path_dataset_train', required=True, type=str)
parser.add_argument('--input_path_dataset_test',  required=True, type=str)
parser.add_argument('--input_path_models',        required=True, type=str)
parser.add_argument('--selected_model',           required=True, type=str)
parser.add_argument('--cv',                       required=True, type=int)
parser.add_argument('--optimizer',                required=True, type=str)
parser.add_argument('--n_iter',                   required=True, type=int)
parser.add_argument('--n_jobs',                   required=True, type=int)
parser.add_argument('--output_dir_models',            required=True, type=str)
parser.add_argument('--output_path_hyperparams',      required=True, type=str)
parser.add_argument('--output_path_performances',     required=True, type=str)
parser.add_argument('--output_path_selected_features',required=True, type=str)
parser.add_argument('--output_path_shap',             required=True, type=str)
parser.add_argument('--tx_id',                    required=True, type=str)
args = parser.parse_args()

CV        = args.cv
OPTIMIZER = args.optimizer
N_ITER    = args.n_iter
N_JOBS    = args.n_jobs
TxID      = args.tx_id
SEED      = 48047510

METRICS_LIST = ["r2", "adj_r2", "explained_variance", "neg_mean_squared_error",
                "neg_root_mean_squared_error", "neg_mean_absolute_error",
                "neg_median_absolute_error", "neg_max_error",
                "d2_absolute_error_score"]

# ---- load pre-split matrices ----
TRAIN = pd.read_csv(args.input_path_dataset_train, header=0, sep="\t", index_col=0)
TEST  = pd.read_csv(args.input_path_dataset_test,  header=0, sep="\t", index_col=0)

assert TxID in TRAIN.columns and TxID in TEST.columns, f"{TxID} missing from a matrix"
assert list(TRAIN.columns) == list(TEST.columns), "train/test columns differ"

feature_cols = [c for c in TRAIN.columns if c != TxID]
X_train, y_train = TRAIN[feature_cols], TRAIN[TxID]
X_test,  y_test  = TEST[feature_cols],  TEST[TxID]

assert not (set(X_train.index) & set(X_test.index)), "train/test sample overlap"

n_feat = X_train.shape[1]
print(f"[{TxID}] train {X_train.shape[0]}x{n_feat} | test {X_test.shape[0]}x{n_feat}")

models_df = pd.read_csv(args.input_path_models, sep="\t")
MODELS = models_df.model.tolist()
if args.selected_model != "all":
    MODELS = [MODELS[MODELS.index(args.selected_model)]]

best_hyperparams_dict, best_models_dict = {}, {}
selected_features_dict, performances_dict, shap_dict = {}, {}, {}

script_start = time.time()

for MODEL in MODELS:
    loop_start = time.time()

    NEEDS_SCALING = bool(models_df[models_df.model == MODEL].needs_scaling.values)
    IS_TREE_BASED = bool(models_df[models_df.model == MODEL].is_tree_based.values)
    REGRESSOR     = BUILD_MODEL(MODEL)

    # k-grid + model params, capped at this transcript's feature count
    HYPERPARAM_SPACE = GET_HYPERPARAM_SPACE_KBEST(MODEL, n_feat)

    print(MODEL)

    t0 = time.time()
    steps = []
    if NEEDS_SCALING:
        steps.append(("scaling", sk.preprocessing.StandardScaler()))
    # cheap univariate selection; k is tuned by the search (leakage-clean: re-fit per fold)
    steps.append(("feature_selection", SelectKBest(score_func=f_regression)))
    steps.append(("model", sk.base.clone(REGRESSOR)))
    pipe = sk.pipeline.Pipeline(steps)
    print(f"[{MODEL}] Pipeline built in {time.time() - t0:.2f}s")

    t0 = time.time()
    best_model, best_hyperparams, cv_score = OPTIMIZE_HYPERPARAMS(
    pipe, HYPERPARAM_SPACE, OPTIMIZER, X_train, y_train,
    SEED, N_JOBS, N_ITER, CV)
    selected_features = X_train.columns[
    best_model.named_steps["feature_selection"].get_support()]
    print(f"[{MODEL}] Hyperparameter optimization done in {time.time() - t0:.2f}s")
    print(f"[{MODEL}] Selected {len(selected_features)} features "
    f"(k={best_hyperparams.get('feature_selection__k')})")

    t0 = time.time()
    performances = {"cv_r2": cv_score}
    performances.update(EVALUATE_PERFORMANCES(best_model, X_train, y_train, METRICS_LIST, prefix="train"))
    performances.update(EVALUATE_PERFORMANCES(best_model, X_test,  y_test,  METRICS_LIST, prefix="test"))
    print(f"[{MODEL}] Performance evaluation done in {time.time() - t0:.2f}s")

    model_final = best_model.named_steps["model"]
    X_train_transformed = best_model[:-1].transform(X_train)

    t0 = time.time()
    explainer = GET_EXPLAINER(model_final, IS_TREE_BASED, X_train_transformed)
    if isinstance(explainer, sh.KernelExplainer):
        background = sh.sample(X_train_transformed, min(50, X_train_transformed.shape[0]))
        explainer = sh.KernelExplainer(model_final.predict, background)
        shap_vals = np.array(explainer.shap_values(X_train_transformed, nsamples=100))
    else:
        shap_vals = np.array(explainer(X_train_transformed).values)
    shap_values_df = pd.DataFrame(shap_vals, index=X_train.index, columns=selected_features)
    print(f"[{MODEL}] SHAP values computed in {time.time() - t0:.2f}s")

    best_hyperparams_dict[MODEL]  = best_hyperparams
    best_models_dict[MODEL]       = best_model
    selected_features_dict[MODEL] = selected_features
    performances_dict[MODEL]      = performances
    shap_dict[MODEL]              = shap_values_df
    print(f"[{MODEL}] Total model time: {time.time() - loop_start:.2f}s")

aggregated = pd.DataFrame.from_dict(performances_dict, orient="index")
best_model_name = aggregated["cv_r2"].idxmax()

with open(args.output_path_hyperparams, "w") as fp:
    for k, v in best_hyperparams_dict[best_model_name].items():
        fp.write(f"{k}\t{v}\n")
with open(args.output_path_performances, "w") as fp:
    for k, v in performances_dict[best_model_name].items():
        fp.write(f"{k}\t{v}\n")
with open(args.output_path_selected_features, "w") as fp:
    for feature in selected_features_dict[best_model_name]:
        fp.write(f"{feature}\n")
shap_dict[best_model_name].to_csv(args.output_path_shap, sep="\t", header=True, doublequote=False)

best_obj = best_models_dict[best_model_name]
jbl.dump(best_obj, f"{args.output_dir_models}/{TxID}_{type(best_obj).__name__}.pkl")

print(f"\n=== Total script time: {time.time() - script_start:.2f}s ===")