rule unpack_topK_NTM:
    output:
        "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.bed"
    input:
        "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_filtered.bed"
    log:
        "logs/unpack_topK_NTM/unpack_topK_NTM_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
    shell:
        """
        grep -Fw {wildcards.TxID} {input} > {output} 2> {log}
        """
