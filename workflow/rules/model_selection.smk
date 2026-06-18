rule model_selection:
    output:
        performance_path       = "results/performance/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{test_size}_{CV}_{TxID}",
        sel_feats_path         = "results/sel_feats/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{test_size}_{CV}_{TxID}",
        hyper_path             = "results/hyper/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{test_size}_{CV}_{TxID}",
        shap_path              = "results/shap/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{test_size}_{CV}_{TxID}"
    input:
        dataset = "data/xgb_clustered/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_xgb_input_clustered.txt",
        models  = "config/MODELS"
    params:
        n_jobs    = 10,
        test_size = 0.3,
        cv        = 5,
        n_iter    = 15,
        optimizer = "RandomSearchCV",
        selected_model = "KNeighborsRegressor",
        output_dir_models = "results/models"
    conda:
        "../envs/py_ML.yml"
    log:
        "logs/model_selection/model_selection_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{test_size}_{CV}_{TxID}.log"
    shell:
        """
        mkdir -p {params.output_dir_models}
        PYTHONPATH=. python workflow/scripts/model_selection_time.py \
        --input_path_dataset={input.dataset} \
        --input_path_models={input.models} \
        --selected_model={params.selected_model} \
        --cv={wildcards.CV} \
        --test_size={wildcards.test_size} \
        --optimizer={params.optimizer} \
        --n_iter={params.n_iter} \
        --n_jobs={params.n_jobs} \
        --output_dir_models={params.output_dir_models} \
        --output_path_hyperparams={output.hyper_path} \
        --output_path_performances={output.performance_path} \
        --output_path_selected_features={output.sel_feats_path} \
        --output_path_shap={output.shap_path} \
        --tx_id={wildcards.TxID} 2> {log}
        """
