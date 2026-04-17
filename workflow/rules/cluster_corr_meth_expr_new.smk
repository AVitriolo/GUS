rule cluster_corr_meth_expr:
    output:
        corr = "results/corr_meth_expr/bed/clustered/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}",
        plot = "results/plots/clustered_corr_meth_expr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}.pdf",
    input:
        "results/corr_meth_expr/bed/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}",
    conda:
        "../envs/r_process_corr_meth_expr.yml" 
    log:
        "logs/cluster_corr_meth_expr/cluster_corr_meth_expr_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/cluster_corr_meth_expr_new.R \
        --input_path_corr={input} \
        --output_path_corr={output.corr} \
        --output_path_plot={output.plot} \
        --TxID={wildcards.TxID} 2> {log}
        """
