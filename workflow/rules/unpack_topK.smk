rule unpack_topK:
    output:
        "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{TxID}.bed"
    input:
        "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.bed"
    log:
        "logs/unpack_topK/unpack_topK_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{TxID}"
    shell:
        """
        grep -Fw {wildcards.TxID} {input} > {output} 2>{log}
        """
