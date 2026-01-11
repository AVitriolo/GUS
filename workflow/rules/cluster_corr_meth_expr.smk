rule cluster_corr_meth_expr:
    output:
        corr_pos = "results/corr_meth_expr/bed/clustered/pos_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}",
        corr_neg = "results/corr_meth_expr/bed/clustered/neg_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}",
        plot_pos = "results/plots/clustered_corr_meth_expr/pos_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}.pdf",
        plot_neg = "results/plots/clustered_corr_meth_expr/neg_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}.pdf"   
    input:
        corr_pos = "results/corr_meth_expr/bed/pos_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}",
        corr_neg = "results/corr_meth_expr/bed/neg_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}"
    conda:
        "../envs/r_process_corr_meth_expr.yml" 
    log:
        "logs/cluster_corr_meth_expr/cluster_corr_meth_expr_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/cluster_corr_meth_expr.R \
        --input_path_corr_pos={input.corr_pos} \
        --input_path_corr_neg={input.corr_neg} \
        --output_path_corr_pos={output.corr_pos} \
        --output_path_corr_neg={output.corr_neg} \
        --output_path_plot_neg={output.plot_neg} \
        --output_path_plot_pos={output.plot_pos} \
        --TxID={wildcards.TxID} 2> {log}
        """