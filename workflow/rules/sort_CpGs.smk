rule sort_CpGs:
    output:
        "resources/CpG_data/{current_date}_{assembly_code}_CpGs_coord_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.sorted.bed"
    input:
        "resources/CpG_data/{current_date}_{assembly_code}_CpGs_coord_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.bed"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
	     bedtools sort -i {input} > {output}
        """
        