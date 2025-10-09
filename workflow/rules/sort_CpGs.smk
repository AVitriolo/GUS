rule sort_CpGs:
    output:
        "resources/CpG_data/{current_date}_{assembly_code}_CpGs_coord_{comparison_type}.sorted.bed"
    input:
        "resources/CpG_data/{current_date}_{assembly_code}_CpGs_coord_{comparison_type}.bed"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
	     sortBed -i {input} > {output}
        """
        