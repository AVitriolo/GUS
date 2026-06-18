rule xgb_input_clustered:
    output:
        "data/xgb_clustered/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_xgb_input_clustered.txt"
    input:
        path_rse     = "data/rse_NTM/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}",
        path_mapping = "data/clustering/clusters_maps/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_clusters_map.txt",
        path_counts  = "resources/RNA/kallisto_counts_NTM_{minCount_expr}_{minSamples_expr}.tsv"
    conda:
        "../envs/r_merge_CpGs.yml"
    log:
        "logs/xgb_input_clustered/xgb_input_clustered_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/xgb_input_clustered.R \
        --input_path_rse={input.path_rse} \
        --input_path_mapping={input.path_mapping} \
        --input_path_counts={input.path_counts} \
        --output_path={output} 2> {log}
        """