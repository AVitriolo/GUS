rule model_selection:
    output:
        performances = "results/performance" + TAG + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}_{TxID}",
        sel_feats    = "results/sel_feats" + TAG + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}_{TxID}",
        hyper        = "results/hyper" + TAG + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}_{TxID}",
        shap         = "results/shap" + TAG + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}_{TxID}"
    input:
        train  = "data/xgb_clustered" + TAG + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_xgb_input_clustered_train.txt",
        test   = "data/xgb_clustered" + TAG + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_xgb_input_clustered_test.txt",
        models = "config/MODELS"
    params:
        selected_model = "SVR",
        optimizer      = "RandomSearchCV",
        n_iter         = config["n_iter_rsearch"][0],
        n_jobs         = config["n_jobs_xgboost"][0],
        n_repeats      = config["n_repeats"][0],
        models_dir     = "results/models" + TAG,
        skip_shap      = "--skip_shap" if int(config.get("skip_shap", 0)) else ""
    threads: int(config["n_jobs_xgboost"][0])   # = --n_jobs of the random search, so --cores is not oversubscribed
    conda:
        "../envs/py_ML.yml"
    priority: 50   # finish a gene (unpack -> corr -> clusters -> matrix -> model) before unpacking more: bounds temp RSEs on disk
    log:
        "logs/model_selection/model_selection_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}_{TxID}.log"
    shell:
        """
        mkdir -p {params.models_dir}
        PYTHONPATH=. python workflow/scripts/model_selection_split_kbest_repeated_kcapped.py \
        --input_path_dataset_train={input.train} \
        --input_path_dataset_test={input.test} \
        --input_path_models={input.models} \
        --selected_model={params.selected_model} \
        --cv={wildcards.CV} \
        --n_repeats={params.n_repeats} \
        --optimizer={params.optimizer} \
        --n_iter={params.n_iter} \
        --n_jobs={params.n_jobs} \
        --output_dir_models={params.models_dir} \
        --output_path_hyperparams={output.hyper} \
        --output_path_performances={output.performances} \
        --output_path_selected_features={output.sel_feats} \
        --tx_id={wildcards.TxID} \
        --output_path_shap={output.shap} {params.skip_shap} 2> {log}
        """