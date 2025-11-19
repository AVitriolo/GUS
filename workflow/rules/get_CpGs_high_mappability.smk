rule get_CpGs_high_mappability:
    output:
        "resources/CpG_data/coord/{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.ids"
    input:
        "resources/CpG_data/coord/{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.filtered.bed"
    log:
        "logs/get_CpGs_high_mappability/get_CpGs_high_mappability_{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.log"
    shell:
        """
        cut -f 4 {input} > {output} 2> {log}
        """