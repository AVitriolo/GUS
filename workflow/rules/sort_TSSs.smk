rule sort_TSSs:
    output:
        "resources/reference_data/{current_date}_{assembly_code}_gencode_v{gencode_version}_tss_{tss_subset}.sorted.bed"
    input:
        "resources/reference_data/{current_date}_{assembly_code}_gencode_v{gencode_version}_tss_{tss_subset}.bed"
    conda:
        "/home/Prostate_GAS/workflow/envs/bedtools.yml"
    shell:
        """
	     bedtools sort -i {input} > {output}
        """
        