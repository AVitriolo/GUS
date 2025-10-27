rule train_models:
    output:
        plot_path="results/plots/ML/{current_date}_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{comparison_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_xgb_output_plot_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}.pdf",
        r2_path="results/performance/{current_date}_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{comparison_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_xgb_output_r2_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}",
        obs_pred_path="results/obs_pred/{current_date}_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{comparison_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_xgb_output_obs_pred_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}",
        sel_feats_path="results/sel_feats/{current_date}_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{comparison_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_xgb_output_sel_feats_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}"
    input:
        "data/{current_date}_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{comparison_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_xgbInput_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}"
    conda:
        "../envs/py_ML.yml"
    shell:
        """
	     python workflow/scripts/train_models.py \
        --n_jobs_xgboost=5 \
        --n_jobs_sklearn=1 \
        --test_size=0.3 \
        --cv=10 \
        --n_iter_rsearch=27 \
        --verbosity=2 \
        --num_features_threshold=50 \
        --input_K=100 \
        --hypertune_random_state_rsearch=36 \
        --error_score="raise" \
        --tree_method="hist" \
        --device="cpu" \
        --TxID={wildcards.TxID} \
        --input_path={input} \
        --output_path_plot={output.plot_path} \
        --output_path_r2={output.r2_path} \
        --output_path_obs_pred={output.obs_pred_path} \
        --output_path_sel_feats={output.sel_feats_path}
        """