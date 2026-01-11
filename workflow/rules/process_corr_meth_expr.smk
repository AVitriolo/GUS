rule process_corr_meth_expr:
    output:
        pos = "results/corr_meth_expr/bed/pos_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}",
        neg = "results/corr_meth_expr/bed/neg_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}"
    input:
        corr = "results/corr_meth_expr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{TxID}",
        sel_feats = "results/sel_feats/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}",
        coord = "resources/composite/{assembly_code}_{sample_type}.bed"
    conda: 
        "../envs/r_process_corr_meth_expr.yml"
    log: 
        "logs/process_corr_meth_expr/process_corr_meth_expr_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/process_corr_meth_expr.R \
        --input_path_corr={input.corr} \
        --input_path_sel_feat={input.sel_feats} \
        --input_path_coord={input.coord} \
        --output_path_cor_neg={output.neg} \
        --output_path_cor_pos={output.pos} 2> {log}
        """