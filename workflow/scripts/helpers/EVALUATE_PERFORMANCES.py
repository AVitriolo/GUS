import numpy as np
import sklearn as sk

def adj_r2(REGRESSOR, X, y):
    print("getting adj_r2")
    r2 = sk.metrics.r2_score(y_true =  y, y_pred = REGRESSOR.predict(X))
    n = X.shape[0]
    p = X.shape[1]
    print("r2", r2)
    print(1-(1-r2)*(n-1)/(n-p-1))
    return 1-(1-r2)*(n-1)/(n-p-1)

def EVALUATE_PERFORMANCES(REGRESSOR, X, y, CV, metrics_list):
    print("evaluating performances")
    performances = {}

    for metric in metrics_list:

        if (metric == "adj_r2"):
            evaluator = adj_r2
        else:
            evaluator = sk.metrics.get_scorer(metric)
            print()
        if CV != 1:
            performances[f"train_{metric}"] = sk.model_selection.cross_val_score(REGRESSOR, X, y, scoring=evaluator, cv=CV)
        else:
            performances[f"test_{metric}"] = evaluator(REGRESSOR, X, y)
    print(performances)
    print("performances evaluated")
    return(performances)
