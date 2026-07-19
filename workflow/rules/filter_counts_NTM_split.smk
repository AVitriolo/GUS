rule filter_counts_NTM:
    output:
        counts_train = "resources/RNA/kallisto_counts_NTM_train_{minCount_expr}_{minSamples_expr}.tsv",
        counts_test  = "resources/RNA/kallisto_counts_NTM_test_{minCount_expr}_{minSamples_expr}.tsv",
        TxIDs        = "resources/TxIDs/TxIDs_NTM_{minCount_expr}_{minSamples_expr}",
        tmm_fit      = "resources/RNA/tmm_fit_{minCount_expr}_{minSamples_expr}.rds"
    input:
        NT    = config["RNA_counts_locations_NT"][0],
        M     = config["RNA_counts_locations_M"][0],
        train = "resources/split/train_samples.txt",
        test  = "resources/split/test_samples.txt"
    conda:
        "../envs/r_get_transcriptome.yml"
    log:
        "logs/filter_counts_NTM/filter_counts_NTM_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        Rscript workflow/scripts/filter_counts_NTM_train_test.R \
        --input_expr_values_path_NT={input.NT} \
        --input_expr_values_path_M={input.M} \
        --train_samples_path={input.train} \
        --test_samples_path={input.test} \
        --minCount_expr={wildcards.minCount_expr} \
        --minSamples_expr={wildcards.minSamples_expr} \
        --output_path_counts_train={output.counts_train} \
        --output_path_counts_test={output.counts_test} \
        --output_path_TxIDs={output.TxIDs} \
        --output_path_tmm_fit={output.tmm_fit} 2> {log}
        """