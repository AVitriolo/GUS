rule subset_beta_matrix:
    output:
        temp("data/beta/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_{comparison_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{TxID}")
    input:
        patterns="data/CpGs_by_TxID/{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_{TxID}",
        betas="resources/CpG_data/beta/{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.filtered.txt"
    log:
        "logs/subset_beta_matrix/subset_beta_matrix_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_{comparison_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{TxID}.log"
    shell:
        """
        head -n 1 {input.betas} > {output} 2> {log}
        grep -Fw -f {input.patterns} {input.betas} >> {output} 2>> {log}
        """
