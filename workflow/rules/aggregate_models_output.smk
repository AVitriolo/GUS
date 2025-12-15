rule aggregate_models_output:
    output:
        output_path_hyperparms="results/hyper/aggr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.txt",
        output_path_performances="results/performance/aggr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.txt"
    input:
        check="results/checks/train/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}"
    conda:
        "../envs/r_CpGs.yml"
    log:
        "logs/aggregate_models_output/aggregate_models_output_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        Rscript workflow/scripts/aggregate_models_output.R \
        --input_dir_hyperparams='results/hyper/' \
        --output_path_hyperparams={output.output_path_hyperparms} \
        --input_dir_performances='results/performance/' \
        --output_path_performances={output.output_path_performances} 2> {log}
        """
