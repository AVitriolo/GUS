rule xgb_input_clustered_split:
    output:
        train = "data/xgb_clustered/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_xgb_input_clustered_train.txt",
        test  = "data/xgb_clustered/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_xgb_input_clustered_test.txt"
    input:
        rse          = "data/rse_NTM/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}",
        mapping      = "data/clustering/clusters_maps/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_clusters_map.txt",
        counts_train = counts_train,
        counts_test  = counts_test,
        train        = "resources/split/train_samples.txt",
        test         = "resources/split/test_samples.txt"
    conda:
        "../envs/r_merge_CpGs.yml"
    log:
        "logs/xgb_input_clustered/xgb_input_clustered_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/xgb_input_clustered_split.R \
        --input_path_rse={input.rse} \
        --input_path_mapping={input.mapping} \
        --input_path_counts_train={input.counts_train} \
        --input_path_counts_test={input.counts_test} \
        --train_samples={input.train} \
        --test_samples={input.test} \
        --output_path_train={output.train} \
        --output_path_test={output.test} 2> {log}
        """