rule unpack_RSE:
    output:
        directory("data/checks/rse/{assembly_code}_{sample_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{TxID}")
    input:
        input_dir="resources/cpgea_wgbs_with_coverage_{assembly_code}",
        path_filtered_CpGs="resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{TxID}.bed"
    conda:
        "../envs/r_xgb_input.yml"
    log:
        "logs/unpack_RSE/unpack_RSE_{assembly_code}_{sample_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/unpack_RSE.R \
        --input_dir={input.input_dir} \
        --input_path_filtered_CpGs={input.path_filtered_CpGs} \
        --output_path_rse={output} \
        --sample_type={wildcards.sample_type} 2> {log}
        """

