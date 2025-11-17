rule get_transcriptome:
    output:
        temp("data/RNA/{sample_type}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{TxID}")
    input:
        input_filtered_expr_values_path="resources/RNA/kallisto_counts_{minCount_expr}_{minSamples_expr}.tsv",
    conda:
        "../envs/r_get_transcriptome.yml"
    log:
        "logs/get_transcriptome/get_transcriptome_{sample_type}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/get_transcriptome.R \
        input_filtered_expr_values_path={input.input_filtered_expr_values_path} \
        TxID={wildcards.TxID} \
        sample_type={wildcards.sample_type} \
        output_path={output} 2> {log}
        """