rule filter_counts:
    output:
        counts="resources/RNA/kallisto_counts_{minCount_expr}_{minSamples_expr}.tsv",
        TxIDs="resources/TxIDs/TxIDs_{minCount_expr}_{minSamples_expr}"
    input:
        "resources/RNA/kallisto_counts.tsv"
    conda:
        "../envs/r_get_transcriptome.yml"
    log:
        "logs/filter_counts/filter_counts_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        Rscript workflow/scripts/filter_counts_new.R \
        --input_expr_values_path={input} \
        --minCount_expr={wildcards.minCount_expr} \
        --minSamples_expr={wildcards.minSamples_expr} \
        --output_path_counts={output.counts} \
        --output_path_TxIDs={output.TxIDs} 2> {log}
        """