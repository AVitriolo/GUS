rule generate_feature_heatmap:
    output:
        "results/plots/corr/{current_date}_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_{comparison_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}"
    input:
        "data/corr/{current_date}_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_{comparison_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}"
    conda:
        "../envs/r_heatmap.yml"
    shell:
        """
        Rscript workflow/scripts/generate_feature_heatmap.R \
        --input_path={input} \
        --output_path={output}        
        """
