rule generate_feature_heatmap:
    output:
        "results/plots/corr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}.pdf"
    input:
        "data/corr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}"
    conda:
        "../envs/r_heatmap.yml"
    log:
        "logs/generate_feature_heatmap/generate_feature_heatmap_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/generate_feature_heatmap.R \
        --input_path={input} \
        --output_path={output} 2> {log}
        """
