rule model_selection:
    output:
        performances = "results/performance/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}_{TxID}",
        sel_feats    = "results/sel_feats/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}_{TxID}",
        hyper        = "results/hyper/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}_{TxID}",
        shap         = "results/shap/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}_{TxID}"
    input:
        train  = "data/xgb_clustered/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_xgb_input_clustered_train.txt",
        test   = "data/xgb_clustered/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_xgb_input_clustered_test.txt",
        models = "config/MODELS"
    params:
        selected_model = "SVR",
        optimizer      = "RandomSearchCV",
        n_iter         = config["n_iter_rsearch"][0],
        n_jobs         = config["n_jobs_xgboost"][0]
    conda:
        "../envs/py_ML.yml"
    log:
        "logs/model_selection/model_selection_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}_{TxID}.log"
    shell:
        """
        mkdir -p results/models
        PYTHONPATH=. python workflow/scripts/model_selection_split_kbest.py \
        --input_path_dataset_train={input.train} \
        --input_path_dataset_test={input.test} \
        --input_path_models={input.models} \
        --selected_model={params.selected_model} \
        --cv={wildcards.CV} \
        --optimizer={params.optimizer} \
        --n_iter={params.n_iter} \
        --n_jobs={params.n_jobs} \
        --output_dir_models=results/models \
        --output_path_hyperparams={output.hyper} \
        --output_path_performances={output.performances} \
        --output_path_selected_features={output.sel_feats} \
        --tx_id={wildcards.TxID} \
        --output_path_shap={output.shap} 2> {log}
        """