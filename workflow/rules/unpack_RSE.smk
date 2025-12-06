rule unpack_RSE:
    output:
        "results/checks/rse/{assembly_code}_{sample_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}"
    input:
        input_dir="resources/cpgea_wgbs_with_coverage_{assembly_code}",
        output_dir="data/rse",
        ff_selected_TxIDs="resources/TxIDs/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_ff_selected_TxIDs",
        path_filtered_CpGs="resources/CpG_data/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.bed"
    conda:
        "../envs/r_xgb_input.yml"
    log:
        "logs/unpack_RSE/unpack_RSE_{assembly_code}_{sample_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        touch {output}
        Rscript workflow/scripts/unpack_RSE.R \
        --input_dir={input.input_dir} \
        --output_dir={input.output_dir} \
        --input_path_selected_TxIDs={input.ff_selected_TxIDs} \
        --input_path_filtered_CpGs={input.path_filtered_CpGs} \
        --sample_type={wildcards.sample_type} 2> {log}
        """
    