rule get_transcriptome:
    output:
        temp("data/{current_date}_{sample_type}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{TxID}_transcriptome")
    input:
        input_filtered_expr_values_path="resources/RNA/kallisto_counts_filtered_{minCount_expr}_{minSamples_expr}.tsv",
    conda:
        "../envs/r_get_transcriptome.yml"
    shell:
        """
        Rscript workflow/scripts/get_transcriptome.R \
        input_filtered_expr_values_path={input.input_filtered_expr_values_path} \
        TxID={wildcards.TxID} \
        sample_type={wildcards.sample_type} \
        output_path={output}
        """