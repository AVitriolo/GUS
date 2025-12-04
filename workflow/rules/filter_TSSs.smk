rule filter_TSSs:
    output:
        filtered_TSSs="resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}.sorted.filtered.bed"
    input:
        unfiltered_TSSs="resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}.sorted.bed",
        filtered_selected_TxIDs="resources/TxIDs/{minCount_expr}_{minSamples_expr}_f_selected_TxIDs"
    log:
        "logs/filter_TSSs/filter_TSSs_{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        grep -Fw -f {input.filtered_selected_TxIDs} {input.unfiltered_TSSs} > {output.filtered_TSSs} 2> {log}
        """
