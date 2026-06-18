rule compute_pairwise_correlation:
    output:
        "data/corr_NTM/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}"
    input:
        path_rse = "data/rse_NTM/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}"
    conda:
        "../envs/r_xgb_input.yml"
    log:
        "logs/corr_NTM/corr_NTM_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/compute_pairwise_correlation.R \
        --input_path_rse={input.path_rse} \
        --output_path_corr={output} 2> {log}
        """