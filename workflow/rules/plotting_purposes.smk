rule plotting_purposes:
    output:
        path_nTxIDxCpG="results/plots/surviving/nTxIDxCpG/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.pdf",
        path_nCpGxTxID="results/plots/surviving/nCpGxTxID/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.pdf",
        path_distance="results/plots/surviving/distance/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.pdf"
    input:
        "results/surviving/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}"
    log:
        "logs/plotting_purposes/plotting_purposes_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}.log"
    conda:
        "../envs/r_heatmap.yml"
    shell:
        """
        Rscript workflow/scripts/plotting_purposes.R \
        --input_path={input} \
        --output_path_nTxIDxCpG={output.path_nTxIDxCpG} \
        --output_path_nCpGxTxID={output.path_nCpGxTxID} \
        --output_path_distance={output.path_distance} 2>{log}
        """