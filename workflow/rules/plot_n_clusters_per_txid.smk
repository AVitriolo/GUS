rule plot_n_clusters_per_txid:
    output:
        "results/plots/clustering/n_clusters_per_txid/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.pdf"
    input:
        get_txid_merge_CpGs_NTM
    params:
        input_dir = "data/clustering/clusters_maps"
    conda:
        "../envs/r_merge_CpGs.yml"
    log:
        "logs/plot_n_clusters_per_txid/plot_n_clusters_per_txid_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        Rscript workflow/scripts/plot_n_clusters_per_txid.R \
        --input_dir_maps={params.input_dir} \
        --output_path_n_clusters={output} 2> {log}
        """