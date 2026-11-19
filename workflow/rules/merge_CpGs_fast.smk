rule merge_CpGs_fast:   # = compute_pairwise_correlation + merge_CpGs in one step, same cluster map (see the script header)
    output:
        path_map  = "data/clustering/clusters_maps" + TAG + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_clusters_map.txt",
        path_plot = "results/plots/clustering" + TAG + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_heatmap.pdf"
    input:
        path_rse    = "data/rse_NTM/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}",
        path_counts = counts_train,
        train       = "resources/split/train_samples.txt"
    params:
        C = 2,
        j = 1,
        plot = int(config.get("clustering_heatmap", 1))   # 0: placeholder PDF (the split heatmap is most of the runtime)
    threads: 2
    conda:
        "../envs/r_merge_CpGs.yml"
    priority: 30   # finish a gene (unpack -> clusters -> matrix -> model) before unpacking more: bounds temp RSEs on disk
    log:
        "logs/merge_CpGs/merge_CpGs_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
    shell:
        """
        OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript workflow/scripts/merge_CpGs_fast.R \
        --input_path_rse={input.path_rse} \
        --input_path_counts={input.path_counts} \
        --train_samples={input.train} \
        --output_path_plot={output.path_plot} \
        --output_path_map={output.path_map} \
        --tx_id={wildcards.TxID} \
        --C={params.C} \
        --j={params.j} \
        --plot={params.plot} 2> {log}
        """
