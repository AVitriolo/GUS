rule get_CpGs_by_TxID:
    output:
        temp("data/CpGs_by_TxID/{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_{TxID}")
    input:
        "data/coord/{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.distance.filtered.bed"
    log:
        "logs/get_CpGs_by_TxID/get_CpGs_by_TxID_{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_{TxID}.log"
    shell:
        """
        grep -Fw {wildcards.TxID} {input} | cut -f 1 > {output} 2> {log} || echo "No matching CpGs found or {wildcards.TxID} is not in {input} >> {log}"
        """