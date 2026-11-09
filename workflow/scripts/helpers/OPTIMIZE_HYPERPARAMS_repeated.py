import sklearn as sk
import skopt as sko
from sklearn.model_selection import RepeatedKFold

# Ablation of OPTIMIZE_HYPERPARAMS.py: pass a RepeatedKFold splitter instead of a
# bare fold-count int. sklearn resolves a plain int `cv` on a regressor to a
# single, non-repeated KFold - fine when you have plenty of samples, but at
# ~239 train samples / 3 folds (~80/fold) that single split is itself a real
# source of variance in cv_r2, on top of R2's unboundedness below on small
# folds. Repeating the split N_REPEATS times and averaging (what
# best_score_ does automatically over all CV*N_REPEATS folds) attacks that
# specific source of instability without needing more data.
def OPTIMIZE_HYPERPARAMS_repeated(ESTIMATOR, HYPERPARAM_SPACE, OPTIMIZER, X, y, SEED, N_JOBS, N_ITER, CV, N_REPEATS):
    print("starting hyperparams optimization")
    HYPERPARAM_SPACE = HYPERPARAM_SPACE[OPTIMIZER]

    cv_splitter = RepeatedKFold(n_splits=CV, n_repeats=N_REPEATS, random_state=SEED)
    print(f"cv splitter: RepeatedKFold(n_splits={CV}, n_repeats={N_REPEATS}) "
          f"-> {CV * N_REPEATS} total fold evaluations per candidate")

    # fix space for BayesSearch
    if OPTIMIZER == "GridSearchCV":

        opti = sk.model_selection.GridSearchCV(estimator = ESTIMATOR,
                                               param_grid = HYPERPARAM_SPACE,
                                               scoring = "r2",
                                               cv = cv_splitter,
                                               refit = True,
                                               n_jobs = N_JOBS,
                                               verbose = 0)

    elif OPTIMIZER == "RandomSearchCV":

        opti = sk.model_selection.RandomizedSearchCV(estimator = ESTIMATOR,
                                                     param_distributions = HYPERPARAM_SPACE,
                                                     n_iter = N_ITER,
                                                     scoring = "r2",
                                                     cv = cv_splitter,
                                                     random_state = SEED,
                                                     refit = True,
                                                     n_jobs = N_JOBS,
                                                     verbose = 0)

    elif OPTIMIZER == "BayesSearchCV":

        opti = sko.BayesSearchCV(estimator = ESTIMATOR,
                                 search_spaces = HYPERPARAM_SPACE,
                                 n_iter = N_ITER,
                                 scoring = "r2",
                                 cv = cv_splitter,
                                 random_state = SEED,
                                 refit = True,
                                 n_jobs = N_JOBS,
                                 verbose = 0)

    else:

        pass

    opti.fit(X, y)
    best_model = opti.best_estimator_
    best_hyperparams = opti.best_params_
    best_cv_score = opti.best_score_  # mean R2 across CV*N_REPEATS folds, used for model selection
    print(best_hyperparams)
    print(f"best CV r2: {best_cv_score}")
    print("done hyperparams optimization")
    return([best_model, best_hyperparams, best_cv_score])
