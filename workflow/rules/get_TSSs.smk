rule get_TSSs:
    output:
	    "resources/reference_data/{current_date}_{assembly_code}_gencode_v{gencode_version}_tss_{tss_subset}.bed"
    input:
        "resources/reference_data/gencode.v{gencode_version}.annotation.gff3.gz"
    conda:
        "../envs/r_TSSs.yml"
    shell:
        """
         Rscript workflow/scripts/get_TSSs.R \
        --input_path={input} \
        --output_path={output}
        """