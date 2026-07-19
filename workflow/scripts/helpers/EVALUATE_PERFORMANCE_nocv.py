import numpy as np
import sklearn as sk

def adj_r2(REGRESSOR, X, y):
    r2 = REGRESSOR.score(X, y)
    n = X.shape[0]
    # count features the model actually uses, not the raw input width
    if hasattr(REGRESSOR, "named_steps") and "feature_selection" in REGRESSOR.named_steps:
        p = int(REGRESSOR.named_steps["feature_selection"].get_support().sum())
    else:
        p = X.shape[1]
    denom = n - p - 1
    return 1 - (1 - r2) * (n - 1) / denom

def EVALUATE_PERFORMANCES(REGRESSOR, X, y, metrics_list, prefix="train"):
    performances = {}

    for metric in metrics_list:

        if metric == "adj_r2":
            evaluator = adj_r2
        else:
            evaluator = sk.metrics.get_scorer(metric)

        performances[f"{prefix}_{metric}"] = evaluator(REGRESSOR, X, y)

    return performances
