import numpy as np
import pandas as pd
import math as mh
import scipy as sp
import sklearn as sk
import skopt as sko
import xgboost as xgb
import shap as sh
import pickle as pkl
import joblib as jbl
import seaborn as sns
import matplotlib.pyplot as plt
import random
import warnings
import argparse

from workflow.scripts.helpers.BUILD_MODEL import *
from workflow.scripts.helpers.GET_HYPERPARAM_SPACE import *

from workflow.scripts.helpers.SPLIT_DATASET import *
from workflow.scripts.helpers.CHECK_CHR_SORTING import *
from workflow.scripts.helpers.BLOCK_DATASET import *

from workflow.scripts.helpers.GENERATE_KFOLD_IDXS import *
from workflow.scripts.helpers.SUBSET_DATA_BY_IDXS import *
from workflow.scripts.helpers.OPTIMIZE_HYPERPARAMS import *
from workflow.scripts.helpers.GET_EXPLAINER import *

from workflow.scripts.helpers.EVALUATE_PERFORMANCES import *
from workflow.scripts.helpers.AGGREGATE_SELECTED_FEATURES import *

parser = argparse.ArgumentParser()

parser.add_argument('--input_path_dataset', dest='input_path_dataset', required=True, type=str, help='')
parser.add_argument('--input_path_models', dest='input_path_models', required=True, type=str, help='')
parser.add_argument('--selected_model', dest = 'selected_model', required = True, type = str, help = '')
parser.add_argument('--cv', dest = 'cv', required = True, type = int, help = '')
parser.add_argument('--test_size', dest = 'test_size', required = True, type = float, help = '')
parser.add_argument('--optimizer', dest = 'optimizer', required = True, type = str, help = '')
parser.add_argument('--n_iter', dest = 'n_iter', required = True, type = int, help = '')
parser.add_argument('--n_jobs', dest = 'n_jobs', required = True, type = int, help = '')
parser.add_argument('--output_dir_models', dest = 'output_dir_models', required = True, type = str, help = '')
parser.add_argument('--output_path_hyperparams', dest = 'output_path_hyperparams', required = True, type = str, help = '')
parser.add_argument('--output_path_performances', dest = 'output_path_performances', required = True, type = str, help = '')
parser.add_argument('--output_path_selected_features', dest = 'output_path_selected_features', required = True, type = str, help = '')
parser.add_argument('--output_path_shap', dest = 'output_path_shap', required = True, type = str, help = '')
parser.add_argument('--tx_id', dest = 'tx_id', required = True, type = str, help = '')

args = parser.parse_args()

DATASET = pd.read_csv(args.input_path_dataset, header = 0, sep = "\t", index_col=0) #!# we iterate over datasets
TEST_SIZE = args.test_size
CV = args.cv
OPTIMIZER = args.optimizer
N_ITER = args.n_iter
N_JOBS = args.n_jobs
TxID = args.tx_id
SEED = 48047510 #!# TODO: list of seeds

METRICS_LIST = ["r2", "adj_r2", "explained_variance", "neg_mean_squared_error",
                "neg_root_mean_squared_error", "neg_mean_absolute_error", "neg_median_absolute_error",
                "neg_max_error", "d2_absolute_error_score"]

models_df = pd.read_csv(args.input_path_models, sep="\t")
MODELS = models_df.model.tolist()

if args.selected_model != "all":
    MODELS = [MODELS[MODELS.index(args.selected_model)]]

# with open(args.input_path_seeds) as f:
#   SEEDS = [line.rstrip() for line in f.readlines()]

num_rows_dataset = DATASET.shape[0]

best_hyperparams_dict = dict()
best_models_dict = dict()
selected_features_dict = dict()
performances_dict = dict()
shap_dict = dict()


for MODEL in MODELS:
    HANDLES_NAs = bool(models_df[models_df.model == MODEL].handles_NAs.values)
    NEEDS_SCALING = bool(models_df[models_df.model == MODEL].needs_scaling.values)
    IS_TREE_BASED = bool(models_df[models_df.model == MODEL].is_tree_based.values)
    REGRESSOR = BUILD_MODEL(MODEL)
    HYPERPARAM_SPACE = GET_HYPERPARAM_SPACE(MODEL) # add hyperparam spaces
    bool_optimize_hyperparams = True

    X_train, X_test, \
    y_train, y_test = SPLIT_DATASET(DATASET, SEED, TEST_SIZE, TxID)

    num_rows_test_set = X_test.shape[0]
    empirical_test_size = (num_rows_test_set / num_rows_dataset)
    
    row_name = f"{MODEL}"
    print(row_name)

    # if model handles NAs --> do not impute and do not scale (because it's XGBRegressor)
    # else --> first impute, then scale 
    if (not NEEDS_SCALING):

        pipe = sk.pipeline.Pipeline([
            ("feature_selection", sk.feature_selection.SequentialFeatureSelector(estimator = sk.base.clone(REGRESSOR), n_features_to_select = "auto", direction = 'forward', cv = CV, n_jobs = N_JOBS)),
            ("model", sk.base.clone(REGRESSOR))
        ])

    else:

        pipe = sk.pipeline.Pipeline([
            ("scaling", sk.preprocessing.StandardScaler()),
            ("feature_selection", sk.feature_selection.SequentialFeatureSelector(estimator = sk.base.clone(REGRESSOR), n_features_to_select = "auto", direction = 'forward', cv = CV, n_jobs = N_JOBS)),
            ("model", sk.base.clone(REGRESSOR))
        ])

    best_model, best_hyperparams = OPTIMIZE_HYPERPARAMS(pipe, HYPERPARAM_SPACE, OPTIMIZER, X_train, y_train, SEED, N_JOBS, N_ITER, CV)
    selected_features = X_train.columns[best_model.named_steps["feature_selection"].get_support()]

    performances_train = EVALUATE_PERFORMANCES(best_model, X_train, y_train, CV, METRICS_LIST)
    performances_test = EVALUATE_PERFORMANCES(best_model, X_test, y_test, 1, METRICS_LIST)

    performances = dict()
    for performance in [performances_train, performances_test]:
         performances.update(performance)
    
    model_final = best_model.named_steps["model"]
    print(model_final)

    X_train_transformed = best_model[:-1].transform(X_train)
    
    explainer = GET_EXPLAINER(model_final, IS_TREE_BASED, X_train_transformed)

    if isinstance(explainer, sh.KernelExplainer):
        shap_vals = np.array(explainer.shap_values(X_train_transformed)) 
    else:
        shap_exp = explainer(X_train_transformed)
        shap_vals = np.array(shap_exp.values)
    
    shap_values_df = pd.DataFrame(shap_vals, index = X_train.index, columns = selected_features)

    best_hyperparams_dict[row_name] = best_hyperparams
    best_models_dict[row_name] = best_model
    selected_features_dict[row_name] = selected_features
    performances_dict[row_name] = performances
    shap_dict[row_name] = shap_values_df
    print("done ", MODEL)

aggregated_performances = pd.DataFrame.from_dict(performances_dict, orient="index")
best_model_name = aggregated_performances["r2"].idxmax()

best_model_best_hyperparameters = best_hyperparams_dict[best_model_name]
best_model_selected_features = selected_features_dict[best_model_name]
best_model_performances = performances_dict[best_model_name]
best_model_object = best_models_dict[best_model_name]
best_model_shap = shap_dict[best_model_name]

print(best_model_best_hyperparameters)
print(best_model_selected_features)
print(best_model_performances)
print(best_model_object)
print(best_model_shap)

with open(args.output_path_hyperparams, 'w') as fp:
        for k,v in best_model_best_hyperparameters.items():
            row = k + "\t" + str(v) + "\n"
            fp.write(row)

with open(args.output_path_performances, 'w') as fp:
        for k,v in best_model_performances.items():
            row = k + "\t" + str(v) + "\n"
            fp.write(row)

with open(args.output_path_selected_features, 'w') as fp:
        for feature in best_model_selected_features:
            fp.write(str(feature) + "\n")

best_model_shap.to_csv(args.output_path_shap, sep = "\t", header=True, doublequote=False)

jbl.dump(best_model_object, (args.output_dir_models + "/" + type(best_model_object).__name__ + ".pkl"))

