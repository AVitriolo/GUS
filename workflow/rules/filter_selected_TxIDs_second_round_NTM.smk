checkpoint filter_selected_TxIDs_second_round_NTM:
    output:
        ff_selected_TxIDs  = "resources/TxIDs/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_ff_selected_TxIDs",
        path_filtered_CpGs = "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_filtered.bed",
        path_num_CpGs_plot = "results/plots/filtering/nCpGs/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.pdf",
        path_dist_filt_plot= "results/plots/filtering/distance/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.pdf"
    input:
        path_CpGs = "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.bed"
    conda:
        "../envs/r_xgb_input.yml"
    log:
        "logs/filter_selected_TxIDs_second_round/filter_selected_TxIDs_second_round_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        Rscript workflow/scripts/filter_selected_TxIDs_second_round_NTM.R \
        --input_path_CpGs={input.path_CpGs} \
        --output_path_ff_selected_TxIDs={output.ff_selected_TxIDs} \
        --output_path_filtered_CpGs={output.path_filtered_CpGs} \
        --output_path_num_CpGs_plot={output.path_num_CpGs_plot} \
        --output_path_dist_filt_plot={output.path_dist_filt_plot} 2> {log}
        """
