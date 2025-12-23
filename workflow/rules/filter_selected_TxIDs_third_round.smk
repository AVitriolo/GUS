checkpoint filter_selected_TxIDs_third_round:
    output:
        fff_selected_TxIDs="resources/TxIDs/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_fff_selected_TxIDs",
    input:
        ff_selected_TxIDs="resources/TxIDs/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_ff_selected_TxIDs",
        TxIDs_to_exclude="resources/TxIDs/TxIDs_to_exclude",
        checks="data/checks/xgb/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}"
    log:
        "logs/filter_selected_TxIDs_third_round/filter_selected_TxIDs_third_round_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}.log"
    shell:
        """
        grep -F -f {input.TxIDs_to_exclude} -v {input.ff_selected_TxIDs} > {output.fff_selected_TxIDs}
        """
