rule unpack_topK_NTM:   # the TxID's candidate CpGs (TxID\tCpGID), pre-split by filter_selected_TxIDs_second_round_NTM
    output:
        "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.bed"
    input:
        "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_split"
    priority: 5   # finish a gene (unpack -> corr -> clusters -> matrix -> model) before unpacking more: bounds temp RSEs on disk
    log:
        "logs/unpack_topK_NTM/unpack_topK_NTM_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
    shell:
        """
        cp {input}/{wildcards.TxID}.bed {output} 2> {log}
        """
