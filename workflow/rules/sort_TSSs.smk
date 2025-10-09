rule sort_TSSs:
    output:
        "resources/reference_data/{current_date}_{assembly_code}_gencode_v{gencode_version}_tss_{tss_subset}.sorted.bed"
    input:
        "resources/reference_data/{current_date}_{assembly_code}_gencode_v{gencode_version}_tss_{tss_subset}.bed"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
	     sortBed -i {input} > {output}
        """
        