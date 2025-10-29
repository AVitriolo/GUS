from helpers.do_feature_selection import do_feature_selection
import os
import argparse
import xgboost
import sklearn
import pandas
import numpy
import scipy.stats
import matplotlib.pyplot
import shap

parser = argparse.ArgumentParser()
parser.add_argument('--n_jobs_xgboost', dest='n_jobs_xgboost', type=int, help='Add n_jobs_xgboost')
parser.add_argument('--n_jobs_sklearn', dest='n_jobs_sklearn', type=int, help='Add n_jobs_sklearn')
parser.add_argument('--test_size', dest='test_size', type=float, help='Add test_size')
parser.add_argument('--cv', dest='cv', type=int, help='Add cv')
parser.add_argument('--n_iter_rsearch', dest='n_iter_rsearch', type=int, help='Add n_iter_rsearch')
parser.add_argument('--verbosity', dest='verbosity', type=int, help='Add verbosity')
parser.add_argument('--num_features_threshold', dest='num_features_threshold', type=int, help='Add num_features_threshold')
parser.add_argument('--input_K', dest='input_K', type=int, help='Add input_K')
parser.add_argument('--hypertune_random_state_rsearch', dest='hypertune_random_state_rsearch', type=int, help='Add hypertune_random_state_rsearch')
parser.add_argument('--error_score', dest='error_score', type=str, help='Add error_score')
parser.add_argument('--tree_method', dest='tree_method', type=str, help='Add tree_method')
parser.add_argument('--device', dest='device', type=str, help='Add device')
parser.add_argument('--TxID', dest='TxID', type=str, help='Add TxID')
parser.add_argument('--input_path', dest='input_path', type=str, help='Add input_path')
parser.add_argument('--output_path_plot', dest='output_path_plot', type=str, help='Add output_path_plot')
parser.add_argument('--output_path_r2', dest='output_path_r2', type=str, help='Add output_path_r2')
parser.add_argument('--output_path_obs_pred', dest='output_path_obs_pred', type=str, help='Add output_path_obs_pred')
parser.add_argument('--output_path_sel_feats', dest='output_path_sel_feats', type=str, help='Add output_path_sel_feats')
parser.add_argument('--output_path_feat_imp', dest='output_path_feat_imp', type=str, help='Add output_path_feat_imp')
parser.add_argument('--output_path_hyper', dest='output_path_hyper', type=str, help='Add output_path_hyper')

args = parser.parse_args()

search_space = {
	'learning_rate': scipy.stats.uniform(0.01, 0.29),
	'max_depth': scipy.stats.randint(1, 15),
	'min_child_weight': scipy.stats.randint(1, 51),
	'gamma': scipy.stats.loguniform(0.1, 1.0),
	'subsample': scipy.stats.uniform(0.6, 0.4),
	'colsample_bytree': scipy.stats.uniform(0.5, 0.5),
	'colsample_bylevel': scipy.stats.uniform(0.5, 0.5),
	'reg_alpha': scipy.stats.loguniform(1e-3, 10),
	'reg_lambda': scipy.stats.loguniform(0.1, 100),
	'max_delta_step': scipy.stats.randint(0, 11),
	'n_estimators': scipy.stats.randint(100, 501)
}

xgb_dataset = pandas.read_csv(args.input_path, sep = "\t", header = 0)

bool_feature_selection = xgb_dataset.shape[1] > args.num_features_threshold

random_inits_xgboost = [85] # 76, 4, 51, 8, 39, 72, 65, 83, 21, 85, 34, 20, 86, 47, 96, 63, 77, 39, 3]

for random_init in random_inits_xgboost:

	regressor = xgboost.XGBRegressor(n_jobs = args.n_jobs_xgboost,
			seed = random_init,
			random_state = random_init,
			tree_method = args.tree_method,
			device = args.device)
	
	K = args.input_K if args.input_K < xgb_dataset.shape[0] else xgb_dataset.shape[0] - 2
	X = xgb_dataset.iloc[:,0:K]
	y = xgb_dataset.loc[:,args.TxID]

	random_splits  = [25] # 14, 50, 51, 53, 80, 81, 12, 68, 20, 76, 54, 65, 85, 24, 50, 46, 20, 49, 9]

	for random_split in random_splits:
		
		X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, 
									test_size = args.test_size, 
									random_state = random_split)
		

		optimizer = sklearn.model_selection.RandomizedSearchCV(estimator = regressor, 
					param_distributions = search_space,
					cv = args.cv, 
					n_jobs = args.n_jobs_sklearn,
					n_iter = args.n_iter_rsearch,
					random_state = args.hypertune_random_state_rsearch,
					refit = True,
					verbose = args.verbosity,
					error_score = args.error_score)

		optimizer.fit(X_train, y_train)
		model = optimizer.best_estimator_
		# best_hyperparams = optimizer.best_params_

		# with open(args.output_path_hyper, 'w') as fp:
		# 	for k,v in best_hyperparams.items():
		# 		row = k + "\t" + str(v) + "\n"
		# 		fp.write(row)

		if bool_feature_selection: #do feature selection

			feature_names = model.feature_names_in_
			feature_importance = model.feature_importances_

			output = do_feature_selection(feature_names, 
						feature_importance, 
						model, 
						X_train, 
						X_test, 
						y_train, 
						cv = args.cv,
						output_path_plot=args.output_path_plot,
						plot = True, 
						verbose = True)
			
			if len(output) == 2: # if no feature gets selected, the best prediction XGBoost can make is the average prediction

				avg_observation, selected_features = output[0], output[1]
				predictions = numpy.repeat(avg_observation, y_train.shape[0])

				obs_pred_table = pandas.concat(
					[
				pandas.Series(y_train, name='Observed').reset_index(drop=True),
				pandas.Series(predictions, name='Predicted').reset_index(drop=True)
					], 
					axis=1)

				obs_pred_table.to_csv(args.output_path_obs_pred, sep = "\t", header=True, doublequote=False)

				r2 = sklearn.metrics.r2_score(y_train, predictions)
				with open(args.output_path_r2, 'w') as fp:
					fp.write(str(r2))

				with open(args.output_path_sel_feats, 'w') as fp:
					for feature in selected_features:
						fp.write(str(feature) + "\n")

				with open(args.output_path_feat_imp, 'w') as fp:
    				pass

				exit(0)

			else: # a best subset exists
				X_train_reduced, X_test_reduced, selected_features = output[0], output[1], output[2]
		else: # skip feature selection
			X_train_reduced, X_test_reduced = X_train, X_test
			selected_features = X_train_reduced.columns.tolist()
			matplotlib.pyplot.figure()
			matplotlib.pyplot.text(0.5, 0.5, 'No data available', ha='center', va='center', fontsize=12)
			matplotlib.pyplot.axis('off')
			matplotlib.pyplot.savefig(args.output_path_plot)
			matplotlib.pyplot.close()

		model = model.fit(X_train_reduced, y_train)

		# booster = model.get_booster()
		# importance_weight = booster.get_score(importance_type = "weight")
		# importance_gain = booster.get_score(importance_type = "gain")
		# importance_cover = booster.get_score(importance_type = "cover")

		# explainer = shap.TreeExplainer(model, approximate = False)
		# shap_values = explainer.shap_values(X_train_reduced)
		# features_for_shap = X_train_reduced.columns.tolist()
		# mean_abs_shap = np.mean(np.abs(shap_values), axis = 0).tolist()
		# importance_shap = {features_for_shap[i]:mean_abs_shap[i] for i in range(0, len(features_for_shap))}

		# feature_importance_df = pandas.DataFrame([importance_weight, importance_gain, importance_cover, importance_shap], index = ["weight", "gain", "cover", "shap"])
		# feature_importance_df.to_csv(args.output_path_feat_imp, sep = "\t", header=True, doublequote=False)

		predictions = model.predict(X_train_reduced)

		obs_pred_table = pandas.concat(
			[
		pandas.Series(y_train, name='Observed').reset_index(drop=True),
		pandas.Series(predictions, name='Predicted').reset_index(drop=True)
			], 
			axis=1)

		obs_pred_table.to_csv(args.output_path_obs_pred, sep = "\t", header=True, doublequote=False)

		r2 = sklearn.metrics.r2_score(y_train, predictions)
		with open(args.output_path_r2, 'w') as fp:
			fp.write(str(r2))

		with open(args.output_path_sel_feats, 'w') as fp:
			for feature in selected_features:
				fp.write(str(feature) + "\n")