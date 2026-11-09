import numpy as np
import pandas as pd
import sklearn as sk
import argparse
import time
import sys
import warnings

from sklearn.linear_model import ElasticNet
from sklearn.compose import TransformedTargetRegressor
from sklearn.exceptions import ConvergenceWarning

from sklearn.linear_model import enet_path
from sklearn.model_selection import RepeatedKFold
from joblib import Parallel, delayed
from workflow.scripts.helpers.EVALUATE_PERFORMANCE_nocv import *

# Elastic-net counterpart of model_selection_split_kbest_repeated_kcapped.py.
# Same inputs, same RepeatedKFold (cv x n_repeats, same seed), same metrics and
# output files - only the model changes: StandardScaler -> ElasticNet on ALL
# candidate features, no SelectKBest. The L1 part does the feature selection
# inside the fit (jointly, not one CpG at a time), the L2 part shares weight
# among correlated neighbouring CpGs. The target is standardised too
# (TransformedTargetRegressor) so one alpha grid fits every gene; predictions
# are back-transformed, so R2 is unaffected.
#
# "Selected features" = CpGs with a non-zero coefficient in the refit model.

warnings.filterwarnings("ignore", category=ConvergenceWarning)
sys.stdout = sys.stderr

parser = argparse.ArgumentParser()
parser.add_argument('--input_path_dataset_train', required=True, type=str)
parser.add_argument('--input_path_dataset_test',  required=True, type=str)
parser.add_argument('--cv',                       required=True, type=int)
parser.add_argument('--n_repeats',                required=True, type=int)
parser.add_argument('--n_jobs',                   required=True, type=int)
parser.add_argument('--output_path_hyperparams',      required=True, type=str)
parser.add_argument('--output_path_performances',     required=True, type=str)
parser.add_argument('--output_path_selected_features',required=True, type=str)
parser.add_argument('--tx_id',                    required=True, type=str)
args = parser.parse_args()

TxID = args.tx_id
SEED = 48047510

METRICS_LIST = ["r2", "explained_variance", "neg_mean_squared_error",
                "neg_root_mean_squared_error", "neg_mean_absolute_error",
                "neg_median_absolute_error", "neg_max_error",
                "d2_absolute_error_score"]

TRAIN = pd.read_csv(args.input_path_dataset_train, header=0, sep="\t", index_col=0)
TEST  = pd.read_csv(args.input_path_dataset_test,  header=0, sep="\t", index_col=0)
assert TxID in TRAIN.columns and TxID in TEST.columns, f"{TxID} missing from a matrix"
assert list(TRAIN.columns) == list(TEST.columns), "train/test columns differ"

feature_cols = [c for c in TRAIN.columns if c != TxID]
X_train, y_train = TRAIN[feature_cols], TRAIN[TxID]
X_test,  y_test  = TEST[feature_cols],  TEST[TxID]
assert not (set(X_train.index) & set(X_test.index)), "train/test sample overlap"
print(f"[{TxID}] train {X_train.shape[0]}x{X_train.shape[1]} | test {X_test.shape[0]}x{X_test.shape[1]}")

script_start = time.time()
# Hyperparameter search: same RepeatedKFold splits as the SVR runs, scored by
# mean R2 over all folds. Instead of fitting every (alpha, l1_ratio) pair from
# scratch (GridSearchCV; >10 min/gene at ~50k CpGs), each fold computes the
# whole alpha path per l1_ratio with warm starts (enet_path, what ElasticNetCV
# does internally). Scaling is re-fit inside each fold, as in the pipeline.
L1_RATIOS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
ALPHAS = np.logspace(1, -2.5, 15)  # descending, for warm starts

def fold_scores(tr, va):
    xs = sk.preprocessing.StandardScaler().fit(X_train.values[tr])
    Xtr, Xva = xs.transform(X_train.values[tr]), xs.transform(X_train.values[va])
    mu, sd = y_train.values[tr].mean(), y_train.values[tr].std()
    ytr = (y_train.values[tr] - mu) / sd
    out = np.empty((len(L1_RATIOS), len(ALPHAS)))
    for i, l1 in enumerate(L1_RATIOS):
        _, coefs, _ = enet_path(Xtr, ytr, l1_ratio=l1, alphas=ALPHAS, max_iter=5000)
        pred = Xva @ coefs * sd + mu  # intercept: centred X and y
        out[i] = [sk.metrics.r2_score(y_train.values[va], pred[:, j]) for j in range(len(ALPHAS))]
    return out

splits = list(RepeatedKFold(n_splits=args.cv, n_repeats=args.n_repeats, random_state=SEED).split(X_train))
scores = np.mean(Parallel(n_jobs=args.n_jobs)(delayed(fold_scores)(tr, va) for tr, va in splits), axis=0)
i, j = np.unravel_index(np.argmax(scores), scores.shape)
cv_score = scores[i, j]
best_hyperparams = {"model__regressor__l1_ratio": L1_RATIOS[i], "model__regressor__alpha": ALPHAS[j]}
print(f"best {best_hyperparams}, CV r2 {cv_score:.4f} (search {time.time() - script_start:.0f}s)")

pipe = sk.pipeline.Pipeline([
    ("scaling", sk.preprocessing.StandardScaler()),
    ("model", TransformedTargetRegressor(
        regressor=ElasticNet(alpha=ALPHAS[j], l1_ratio=L1_RATIOS[i], max_iter=5000),
        transformer=sk.preprocessing.StandardScaler())),
])
best_model = pipe.fit(X_train, y_train)

coef = best_model.named_steps["model"].regressor_.coef_
selected_features = X_train.columns[coef != 0]
print(f"Selected {len(selected_features)} features (non-zero coefficients)")

performances = {"cv_r2": cv_score, "n_nonzero": len(selected_features)}
performances.update(EVALUATE_PERFORMANCES(best_model, X_train, y_train, METRICS_LIST, prefix="train"))
performances.update(EVALUATE_PERFORMANCES(best_model, X_test,  y_test,  METRICS_LIST, prefix="test"))

with open(args.output_path_hyperparams, "w") as fp:
    for k, v in best_hyperparams.items():
        fp.write(f"{k}\t{v}\n")
with open(args.output_path_performances, "w") as fp:
    for k, v in performances.items():
        fp.write(f"{k}\t{v}\n")
with open(args.output_path_selected_features, "w") as fp:
    for f, c in zip(selected_features, coef[coef != 0]):
        fp.write(f"{f}\t{c}\n")

print(f"\n=== Total script time: {time.time() - script_start:.2f}s ===")
