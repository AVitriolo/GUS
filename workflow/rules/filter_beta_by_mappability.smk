rule filter_beta_by_mappability:
    output:
        "resources/CpG_data/beta/{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.filtered.txt"
    input:
        IDs="resources/CpG_data/coord/{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.ids",
        beta="resources/CpG_data/beta/{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.txt"
    log:
        "logs/filter_beta_by_mappability/filter_beta_by_mappability_{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.log"
    shell:
        """
        head -n 1 {input.beta} > {output} 2> {log}
        grep -Fw -f {input.IDs} {input.beta} >> {output} 2>> {log}
        """