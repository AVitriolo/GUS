rule filter_beta_cov_plotting_purposes:
    output:
        "results/surviving/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}"
    input:
        path_CpG_mat = config["CpGs_locations"][0],
        path_filtered_top_K = "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.bed"
    log:
        "logs/filter_beta_cov_plotting_purposes/filter_beta_cov_plotting_purposes_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.log"
    conda:
        "../envs/r_xgb_input.yml"
    shell:
        """
        Rscript workflow/scripts/filter_beta_cov_plotting_purposes.R \
        --input_path_CpG_mat={input.path_CpG_mat} \
        --input_path_filtered_top_K={input.path_filtered_top_K} \
        --sample_type={wildcards.sample_type} \
        --minCov={wildcards.minCov} \
        --leftCount_beta={wildcards.leftCount_beta} \
        --rightCount_beta={wildcards.rightCount_beta} \
        --minSamples_beta={wildcards.minSamples_beta} \
        --output_path={output} 2> {log}
        """