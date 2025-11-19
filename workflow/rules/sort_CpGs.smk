rule sort_CpGs:
    output:
        "resources/CpG_data/coord/{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.filtered.sorted.bed"
    input:
        "resources/CpG_data/coord/{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.filtered.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/sort_CpGs/sort_CpGs_{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.log"
    shell:
        """
	     bedtools sort -i {input} > {output} 2> {log}
        """
        