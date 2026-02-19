rule filter_TxIDs_and_CpGs_plotting_purposes:
    output:
        path_filtered_CpGs="resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.bed",
    input:
        path_K_closest="resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}.bed"
    conda:
        "../envs/r_xgb_input.yml"
    log:
        "logs/filter_TxIDs_and_CpGs_plotting_purposes/filter_TxIDs_and_CpGs_plotting_purposes_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        Rscript workflow/scripts/filter_TxIDs_and_CpGs_plotting_purposes.R \
        --input_path_K_closest={input.path_K_closest} \
        --output_path_filtered_CpGs={output.path_filtered_CpGs} \
        --distance={wildcards.distance} \
        --min_CpG={wildcards.min_CpG} 2> {log}
        """