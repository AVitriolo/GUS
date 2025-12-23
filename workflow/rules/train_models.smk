rule train_models:
    output:
        plot_path="results/plots/ML/sel_feats/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}.pdf",
        performance_path="results/performance/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}",
        obs_pred_path="results/obs_pred/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}",
        sel_feats_path="results/sel_feats/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}",
        feat_imp_path="results/feat_imp/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}",
        hyper_path="results/hyper/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}",
        shap_path="results/shap/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}"
    input:
        "data/xgb/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{TxID}"
    params:
        n_jobs_xgboost=config["n_jobs_xgboost"][0],
        n_jobs_sklearn=config["n_jobs_sklearn"][0],
        n_iter_rsearch=config["n_iter_rsearch"][0],
        verbosity=config["verbosity"][0],
        num_features_threshold=config["num_features_threshold"][0],
        hypertune_random_state_rsearch=config["hypertune_random_state_rsearch"][0]
    conda:
        "../envs/py_ML.yml"
    log:
        "logs/train_models/train_models_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}.log"
    shell:
        """
        python workflow/scripts/train_models.py \
        --n_jobs_xgboost={params.n_jobs_xgboost} \
        --n_jobs_sklearn={params.n_jobs_sklearn} \
        --test_size={wildcards.test_size} \
        --cv={wildcards.CV} \
        --n_iter_rsearch={params.n_iter_rsearch} \
        --verbosity={params.verbosity} \
        --num_features_threshold={params.num_features_threshold}\
        --hypertune_random_state_rsearch={params.hypertune_random_state_rsearch} \
        --error_score="raise" \
        --tree_method="hist" \
        --device="cpu" \
        --TxID={wildcards.TxID} \
        --input_path={input} \
        --output_path_plot={output.plot_path} \
        --output_path_performance={output.performance_path} \
        --output_path_obs_pred={output.obs_pred_path} \
        --output_path_sel_feats={output.sel_feats_path} \
        --output_path_feat_imp={output.feat_imp_path} \
        --output_path_hyper={output.hyper_path} \
        --output_path_shap={output.shap_path} 2> {log}
        """

