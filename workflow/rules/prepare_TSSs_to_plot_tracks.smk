rule prepare_TSSs_to_plot_tracks:
    output:
        "resources/reference_data/gencode/{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.sorted.subset.bed"
    input:
        selected_TxIDs="resources/TxIDs/{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_filtered_selected_TxIDs",
        bed="resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}.sorted.bed"
    log:
        "logs/prepare_TSS_to_plot_tracks/prepare_TSS_to_plot_tracks_{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.log"
    shell:
        """
        grep -F -f {input.selected_TxIDs} {input.bed} > {output} 2> {log}
        """