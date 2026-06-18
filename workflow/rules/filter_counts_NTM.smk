rule filter_counts_NTM:
    output:
        counts = "resources/RNA/kallisto_counts_NTM_{minCount_expr}_{minSamples_expr}.tsv",
        TxIDs  = "resources/TxIDs/TxIDs_NTM_{minCount_expr}_{minSamples_expr}"
    input:
        NT = config["RNA_counts_locations_NT"][0],
        M  = config["RNA_counts_locations_M"][0]
    conda:
        "../envs/r_get_transcriptome.yml"
    log:
        "logs/filter_counts_NTM/filter_counts_NTM_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        Rscript workflow/scripts/filter_counts_NTM.R \
        --input_expr_values_path_NT={input.NT} \
        --input_expr_values_path_M={input.M} \
        --minCount_expr={wildcards.minCount_expr} \
        --minSamples_expr={wildcards.minSamples_expr} \
        --output_path_counts={output.counts} \
        --output_path_TxIDs={output.TxIDs} 2> {log}
        """
