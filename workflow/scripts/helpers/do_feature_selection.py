def do_feature_selection(feature_names, feature_importances, model, X_train, X_test, y_train, cv, output_path_plot, plot = True, verbose = True):
	"""
	This function performes incremental feature selection. First takes the feature importance (any metric),
	and trains p (number of features) models, each time adding one feature starting from the least important one
	"""
	import numpy
	import matplotlib.pyplot
	import sklearn

	thresholds = numpy.sort(feature_importances)
	num_features_list = []
	r2_scores_list = []
	for threshold in thresholds:
		vars_to_keep = numpy.where(feature_importances >= threshold)[0]
		X_train_selected = X_train.iloc[:,vars_to_keep]
		r2_scores = sklearn.model_selection.cross_val_score(estimator = model,X = X_train_selected, y = y_train,scoring = "r2",cv = cv)
		mean_r2_score = numpy.mean(r2_scores)
		num_features_list.append(X_train_selected.shape[1])
		r2_scores_list.append(mean_r2_score)
		if verbose:
			print(f'> threshold={threshold}, features={X_train_selected.shape[1]}, R2={mean_r2_score}')
		if len(vars_to_keep) == 1:
			break
	if plot:
		matplotlib.pyplot.figure(figsize=(8, 6))
		matplotlib.pyplot.plot(num_features_list, r2_scores_list, marker='o')
		matplotlib.pyplot.xlabel('Number of Selected Features')
		matplotlib.pyplot.ylabel('Cross Validated mean R2')
		matplotlib.pyplot.title('R2 vs. Number of Selected Features')
		matplotlib.pyplot.grid(True)
		matplotlib.pyplot.savefig(output_path_plot, dpi=300, bbox_inches='tight')
		matplotlib.pyplot.show()
		matplotlib.pyplot.close()
	optimal_threshold_index = numpy.argmax(r2_scores_list)
	optimal_num_features = num_features_list[optimal_threshold_index]
	if (optimal_threshold_index) == 0:
		optimal_threshold = 0
	else:
		optimal_threshold = thresholds[optimal_threshold_index - 1]
	if verbose:
		print(f"Optimal Threshold: {optimal_threshold:.4f}")
		print(f"Number of Selected Features: {optimal_num_features}")
		print(f"R2 at Optimal Threshold: {r2_scores_list[optimal_threshold_index]:.4f}")
		print("\n")
	selected_features = feature_names[numpy.where(feature_importances > optimal_threshold)]
	selected_features = [str(name) for name in selected_features]
	discarded_features = [str(name) for name in feature_names if str(name) not in selected_features]
	if verbose:
		n = 5
		print("Selected features are: \n")
		for i in range(0,len(selected_features), n):
			print("  ".join(selected_features[i:i+n]))
			print("\n")
		print("Discarded features are: \n")
		for i in range(0,len(discarded_features), n):
			print("  ".join(discarded_features[i:i+n]))

	print(selected_features)
	if len(selected_features) > 0:
		X_train_reduced = X_train[selected_features]
		X_test_reduced = X_test[selected_features]
		return([X_train_reduced, X_test_reduced, selected_features])
	else:
		return([numpy.mean(y_train), [""]])
		
	