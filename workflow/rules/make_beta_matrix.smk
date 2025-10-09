rule make_beta_matrix:
    output:
	    "resources/CpG_data/{current_date}_{assembly_code}_CpGs_beta_{sample_type}.txt"
    input:
        "resources/cpgea_wgbs_with_coverage_{assembly_code}"
    conda:
        "../envs/r_CpGs.yml"
    shell:
        """
         Rscript workflow/scripts/make_beta_matrix.R \
        --input_path={input} \
        --output_path={output} \
        --sample_type={wildcards.sample_type}
        """