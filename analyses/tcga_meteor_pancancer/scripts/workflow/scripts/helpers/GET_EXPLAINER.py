import shap as sh

from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import (LinearRegression,Ridge,Lasso,ElasticNet)
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from xgboost import XGBRegressor

def GET_EXPLAINER(MODEL, IS_TREE_BASED, X_BACKGROUND):
    
    if isinstance(MODEL, XGBRegressor):
        explainer = sh.Explainer(MODEL.predict, X_BACKGROUND)

    elif IS_TREE_BASED:
        explainer = sh.TreeExplainer(MODEL,approximate=False)

    elif isinstance(MODEL, (SVR, KNeighborsRegressor)):
        explainer = sh.KernelExplainer(MODEL.predict,X_BACKGROUND)
    
    elif isinstance(MODEL, (LinearRegression, Ridge, Lasso, ElasticNet)):
        explainer = sh.LinearExplainer(MODEL,X_BACKGROUND)

    else:
        explainer = sh.Explainer(MODEL, approximate=False)

    return(explainer)
