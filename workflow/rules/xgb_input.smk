checkpoint xgb_input:
    output:
        path_num_CpGs_plot="results/plots/filtering/nCpGs/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.pdf",
        path_dist_filt_plot="results/plots/filtering/distance/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.pdf",
        ff_selected_TxIDs="resources/TxIDs/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_ff_selected_TxIDs"
    input:
        dir="resources/cpgea_wgbs_with_coverage_{assembly_code}",
        dir_TMRs=config["TMRs_dir"][0],
        path_K_closest="resources/CpG_data/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}.bed",
        path_counts="resources/RNA/kallisto_counts_{minCount_expr}_{minSamples_expr}.tsv"
    params:
        bool_TMRs = config["use_TMRs"][0]
    conda:
        "../envs/r_xgb_input.yml"
    log:
        "logs/xgb_input/xgb_input_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        if [ {params.bool_TMRs} -eq 1 ]; then
            echo "Integrating TMRs"
            script="workflow/scripts/xgb_input_TMRs.R"
            echo "Running $script"
        else
            echo "Not Integrating TMRs"
            script="workflow/scripts/xgb_input.R"
            echo "Running $script"
        fi  

        Rscript $script \
        --input_dir={input.dir} \
        --input_dir_TMRs={input.dir_TMRs} \
        --input_path_K_closest={input.path_K_closest} \
        --input_path_counts={input.path_counts} \
        --output_path_num_CpGs_plot={output.path_num_CpGs_plot} \
        --output_path_dist_filt_plot={output.path_dist_filt_plot} \
        --output_path_ff_selected_TxIDs={output.ff_selected_TxIDs} \
        --sample_type={wildcards.sample_type} \
        --distance={wildcards.distance} \
        --min_CpG={wildcards.min_CpG} \
        --minCov={wildcards.minCov} \
        --leftCount_beta={wildcards.leftCount_beta} \
        --rightCount_beta={wildcards.rightCount_beta} \
        --minSamples_beta={wildcards.minSamples_beta} 2> {log}
        """

