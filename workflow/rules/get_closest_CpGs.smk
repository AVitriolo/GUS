rule get_closest_CpGs:
    output:
        "data/{current_date}_CpGs_coord_{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}.TSS.distance.bed"
    input:
        CpG = "resources/CpG_data/{current_date}_{assembly_code}_CpGs_coord_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}.sorted.bed",
        TSS = "resources/reference_data/{current_date}_{assembly_code}_gencode_v{gencode_version}_tss_{tss_subset}.sorted.bed"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
	     bedtools closest -a {input.CpG} -b {input.TSS} -d | cut -f 4,8,9 > {output}
        """



