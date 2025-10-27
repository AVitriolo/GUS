rule filter_counts:
    output:
        "resources/RNA/kallisto_counts_filtered_{minCount_expr}_{minSamples_expr}.tsv"
    input:
        "resources/RNA/kallisto_counts.tsv"
    conda:
        "../envs/r_get_transcriptome.yml"
    shell:
        """
        Rscript workflow/scripts/filter_counts.R \
        --input_expr_values_path={input} \
        --minCount_expr={wildcards.minCount_expr} \
        --minSamples_expr={wildcards.minSamples_expr} \
        --output_path={output}
        """