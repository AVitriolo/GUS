rule merge_CpGs:
    output:
        path_map  = "data/clustering/clusters_maps/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_clusters_map.txt",
        path_plot = "results/plots/clustering/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_heatmap.pdf"
    input:
        path_corr   = "data/corr_NTM/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}",
        path_rse    = "data/rse_NTM/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}",
        path_counts = counts_train,                                    # <<< train counts
        train       = "resources/split/train_samples.txt"              # <<< NEW
    params:
        C = 4,
        j = 1
    conda:
        "../envs/r_merge_CpGs.yml"
    log:
        "logs/merge_CpGs/merge_CpGs_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/merge_CpGs_split_new.R \
        --input_path={input.path_corr} \
        --input_path_rse={input.path_rse} \
        --input_path_counts={input.path_counts} \
        --train_samples={input.train} \
        --output_path_plot={output.path_plot} \
        --output_path_map={output.path_map} \
        --C={params.C} \
        --j={params.j} 2> {log}
        """