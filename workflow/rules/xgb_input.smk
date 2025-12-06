rule xgb_input:
    output:
        outpath_xgb="data/xgb/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{TxID}",
        outpath_corr="data/corr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{TxID}"
    input:
        "results/checks/rse/{assembly_code}_{sample_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}",
        dir_TMRs=config["TMRs_dir"][0],
        dir_rse="data/rse",
        path_counts="resources/RNA/kallisto_counts_{minCount_expr}_{minSamples_expr}.tsv"
    params:
        bool_TMRs = config["use_TMRs"][0]
    conda:
        "../envs/r_xgb_input.yml"
    log:
        "logs/xgb_input/xgb_input_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
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
        --input_dir_rse={input.dir_rse} \
        --input_dir_TMRs={input.dir_TMRs} \
        --input_path_counts={input.path_counts} \
        --TxID={wildcards.TxID} \
        --sample_type={wildcards.sample_type} \
        --minCov={wildcards.minCov} \
        --leftCount_beta={wildcards.leftCount_beta} \
        --rightCount_beta={wildcards.rightCount_beta} \
        --minSamples_beta={wildcards.minSamples_beta} \
        --output_path_xgb={output.outpath_xgb} \
        --output_path_corr={output.outpath_corr} 2> {log}
        """

