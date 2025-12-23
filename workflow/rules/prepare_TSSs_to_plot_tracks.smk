rule prepare_TSSs_to_plot_tracks:
    output:
        "resources/reference_data/gencode/TSS_to_plot/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}_{TxID}.bed"
    input:
        bed="resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}.sorted.filtered.bed"
    log:
        "logs/prepare_TSS_to_plot_tracks/prepare_TSS_to_plot_tracks_{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
    shell:
        """
        grep -Fw {wildcards.TxID} {input.bed} > {output} 2> {log}
        """