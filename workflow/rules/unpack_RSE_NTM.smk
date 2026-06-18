rule unpack_RSE_NTM:
    output:
        temp(directory("data/rse_NTM/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}"))
    input:
        input_dir          = "resources/integrated_wgbs_hdf5/integrated_wgbs.h5",
        path_filtered_CpGs = "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.bed"
    conda:
        "../envs/r_xgb_input.yml"
    log:
        "logs/unpack_RSE_NTM/unpack_RSE_NTM_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/unpack_RSE_NTM.R \
        --input_dir={input.input_dir} \
        --input_path_filtered_CpGs={input.path_filtered_CpGs} \
        --output_path_rse={output} 2> {log}
        """
