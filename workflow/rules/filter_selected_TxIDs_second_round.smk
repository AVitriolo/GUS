checkpoint filter_selected_TxIDs_second_round:
    output:
        ff_selected_TxIDs="resources/TxIDs/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_ff_selected_TxIDs"
    input:
        path_K_closest="resources/CpG_data/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}.bed"
    conda:
        "../envs/r_xgb_input.yml"
    log:
        "logs/filter_selected_TxIDs_second_round/filter_selected_TxIDs_second_round_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        Rscript workflow/scripts/filter_selected_TxIDs_second_round.R \
        --input_path_K_closest={input.path_K_closest} \
        --output_path_ff_selected_TxIDs={output.ff_selected_TxIDs} \
        --distance={wildcards.distance} \
        --min_CpG={wildcards.min_CpG} 2> {log}
        """
    