checkpoint unpack_RSE:
    output:
        path_num_CpGs_plot="results/plots/filtering/nCpGs/{assembly_code}_{sample_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.pdf",
        path_dist_filt_plot="results/plots/filtering/distance/{assembly_code}_{sample_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.pdf",
    input:
        input_dir="resources/cpgea_wgbs_with_coverage_{assembly_code}",
        output_dir="data/rse",
        path_K_closest="resources/CpG_data/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}.bed",
        ff_selected_TxIDs="resources/TxIDs/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_ff_selected_TxIDs"
    conda:
        "../envs/r_xgb_input.yml"
    log:
        "logs/unpack_RSE/unpack_RSE_{assembly_code}_{sample_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        Rscript workflow/scripts/unpack_RSE.R \
        --input_dir={input.input_dir} \
        --output_dir={input.output_dir} \
        --input_path_selected_TxIDs={input.ff_selected_TxIDs} \
        --input_path_K_closest={input.path_K_closest} \
        --output_path_num_CpGs_plot={output.path_num_CpGs_plot} \
        --output_path_dist_filt_plot={output.path_dist_filt_plot} \
        --sample_type={wildcards.sample_type} \
        --distance={wildcards.distance} \
        --min_CpG={wildcards.min_CpG} 2> {log}
        """
    