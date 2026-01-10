rule get_sel_feats_coordinates:
    output:
        selected_composite_coords="results/composite/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}"
    input:
        CpG_TxID_pairs="results/pairs/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}",
        composite_coords="resources/composite/{assembly_code}_{sample_type}.bed"
    log:
        "logs/get_sel_feats_coordinates/get_sel_feats_coordinates_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}.log"
    shell:
        """
        cat {input.composite_coords} | grep -Fwf <(cut -f 1 {input.CpG_TxID_pairs}) > {output.selected_composite_coords} 2> {log}
        """