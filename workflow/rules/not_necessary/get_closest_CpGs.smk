rule get_closest_CpGs:
    output:
        "data/coord/{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.distance.bed"
    input:
        CpG = "resources/CpG_data/coord/{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.filtered.sorted.bed",
        TSS = "resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}.sorted.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/get_closest_CpGs/get_closest_CpGs_{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.distance.log"
    shell:
        """
	     bedtools closest -a {input.CpG} -b {input.TSS} -d | cut -f 4,8,9 > {output} 2> {log}
        """



