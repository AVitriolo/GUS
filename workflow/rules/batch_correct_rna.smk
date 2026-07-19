rule batch_correct_rna:
    output:
        counts_train = "resources/RNA/kallisto_counts_NTM_train_combat_{minCount_expr}_{minSamples_expr}.tsv",
        counts_test  = "resources/RNA/kallisto_counts_NTM_test_combat_{minCount_expr}_{minSamples_expr}.tsv",
        params_npz   = "resources/RNA/combat_params_{minCount_expr}_{minSamples_expr}.npz",
        plot_ds      = "results/plots/batch_correct/rna_combat_by_dataset_{minCount_expr}_{minSamples_expr}.pdf",
        plot_cond    = "results/plots/batch_correct/rna_combat_by_condition_{minCount_expr}_{minSamples_expr}.pdf"
    input:
        counts_train = "resources/RNA/kallisto_counts_NTM_train_{minCount_expr}_{minSamples_expr}.tsv",
        counts_test  = "resources/RNA/kallisto_counts_NTM_test_{minCount_expr}_{minSamples_expr}.tsv"
    conda:
        "../envs/py_combat.yml"
    log:
        "logs/batch_correct/rna_combat_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        PYTHONPATH=workflow/scripts/helpers python workflow/scripts/combat_correction_rna.py \
        --counts_train={input.counts_train} \
        --counts_test={input.counts_test} \
        --out_train={output.counts_train} \
        --out_test={output.counts_test} \
        --out_params={output.params_npz} \
        --out_plot_ds={output.plot_ds} \
        --out_plot_cond={output.plot_cond} 2> {log}
        """