rule unpack_topK:
    output:
        "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{TxID}.bed"
    input:
        "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}.bed"
    log:
        "logs/unpack_topK/unpack_topK_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{TxID}"
    shell:
        """
        grep -Fw {wildcards.TxID} {input} > {output} 2>{log}
        """
