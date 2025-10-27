rule make_beta_matrix:
    output:
        beta="resources/CpG_data/{current_date}_{assembly_code}_CpGs_beta_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.txt",
        CpGs="resources/CpG_data/{current_date}_{assembly_code}_CpGs_coord_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.bed"
    input:
        "resources/cpgea_wgbs_with_coverage_{assembly_code}"
    conda:
        "../envs/r_CpGs.yml"
    shell:
        """
         Rscript workflow/scripts/make_beta_matrix.R \
        --input_path={input} \
        --output_path_beta={output.beta} \
        --output_path_CpGs={output.CpGs} \
        --sample_type={wildcards.sample_type} \
        --minCov={wildcards.minCov}
        --leftCount_beta={wildcards.leftCount_beta} \
        --rightCount_beta={wildcards.rightCount_beta} \
        --minSamples_beta={wildcards.minSamples_beta}
        """