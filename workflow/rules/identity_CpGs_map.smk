rule identity_CpGs_map:            # cluster_cpgs: 0 -> every candidate CpG is its own feature (replaces merge_CpGs)
    output:
        path_map = "data/clustering/clusters_maps" + TAG + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}_clusters_map.txt"
    input:
        CpGs = "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.bed"
    log:
        "logs/identity_CpGs_map/identity_CpGs_map_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
    shell:
        """
        awk -F'\\t' '{{print $2"\\t"$2}}' {input.CpGs} > {output.path_map} 2> {log}
        """
