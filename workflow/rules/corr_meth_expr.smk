rule corr_meth_expr:
    output:
        "results/corr_meth_expr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{TxID}"
    input:
        "data/xgb/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{TxID}"
    conda:
        "../envs/r_filter_CpGs_TSSs.yml"
    log:
        "logs/corr_meth_expr/corr_meth_expr_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/corr_meth_expr.R \
        --input_path={input} \
        --TxID={wildcards.TxID} \
        --output_path={output} 
        """
