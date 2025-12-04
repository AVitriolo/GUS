rule filter_selected_TxIDs:
    output:
        "resources/TxIDs/{minCount_expr}_{minSamples_expr}_f_selected_TxIDs"
    input:
        selected_TxIDs=config["selected_txIDs_location"][0],
        kallisto_TxIDs="resources/TxIDs/TxIDs_{minCount_expr}_{minSamples_expr}"
    conda:
        "../envs/r_CpGs.yml"
    log:
        "logs/filter_selected_TxIDs/filter_selected_TxIDs_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        Rscript workflow/scripts/filter_selected_TxIDs_new.R \
        --input_path_selected_TxIDs={input.selected_TxIDs} \
        --input_path_kallisto={input.kallisto_TxIDs} \
        --output_path={output} 2> {log}
        """
