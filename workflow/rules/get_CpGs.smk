rule get_CpGs:
    output:
	    "resources/CpG_data/{current_date}_{assembly_code}_CpGs_coord_{comparison_type}.bed"
    input:
        "resources/cpgea_wgbs_with_coverage_{assembly_code}"
    conda:
        "../envs/r_CpGs.yml"
    shell:
        """
         Rscript workflow/scripts/get_CpGs.R \
        --input_path={input} \
        --output_path={output}
        """