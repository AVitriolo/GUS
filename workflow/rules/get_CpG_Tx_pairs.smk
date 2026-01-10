rule get_CpG_Tx_pairs:
    output:
        path_pairs="results/pairs/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}",
    input:
        check="results/checks/train/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}"
    conda:
        "../envs/py_ML.yml"
    log:
        "logs/get_CpG_Tx_pairs/get_CpG_Tx_pairs_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}.log"
    shell:
        """
        python workflow/scripts/get_CpG_Tx_pairs.py \
        --input_dir_selected_CpGs="results/sel_feats" \
        --output_path_pairs={output.path_pairs} 2> {log}
        """