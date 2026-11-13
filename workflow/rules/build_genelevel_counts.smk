checkpoint build_genelevel_counts:  # level: gene -> gene counts (rows = representative TxID) + the representative TxIDs
    output:
        NT             = "resources/RNA/genelevel/kallisto_counts_genelevel_NT_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.tsv",
        M              = "resources/RNA/genelevel/kallisto_counts_genelevel_M_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.tsv",
        representative = "resources/RNA/genelevel/representative_TxIDs_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.tsv"
    input:
        NT               = config["RNA_counts_locations_NT"][0],
        M                = config["RNA_counts_locations_M"][0],
        conversion_table = config["conversion_table"],
        selected_TxIDs   = "resources/TxIDs/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_ff_selected_TxIDs",
        counts_train_tx  = lambda wc: f"resources/RNA/kallisto_counts_NTM_train{'_combat' if BC_RNA else ''}_{wc.minCount_expr}_{wc.minSamples_expr}.tsv"
    conda:
        "../envs/py_ML.yml"
    log:
        "logs/build_genelevel_counts/build_genelevel_counts_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        python workflow/scripts/build_genelevel_counts.py \
        --counts_NT={input.NT} \
        --counts_M={input.M} \
        --conversion_table={input.conversion_table} \
        --selected_TxIDs={input.selected_TxIDs} \
        --counts_train_tx={input.counts_train_tx} \
        --out_NT={output.NT} \
        --out_M={output.M} \
        --out_representative={output.representative} 2> {log}
        """

def genelevel_raw(kind):
    # the gene counts of the (single) configured window / TSS set; the RNA filter rules have no window wildcards
    return (f"resources/RNA/genelevel/kallisto_counts_genelevel_{kind}_{assembly_codes[0]}_v{gencode_versions[0]}_"
            f"{tss_subsets[0]}_{distances[0]}_{{minCount_expr}}_{{minSamples_expr}}.tsv")

rule filter_counts_NTM_genelevel:   # filter_counts_NTM on the gene counts (same script, same 25_10 / TMM)
    output:
        counts_train = "resources/RNA/genelevel/genelevel_counts_NTM_train_{minCount_expr}_{minSamples_expr}.tsv",
        counts_test  = "resources/RNA/genelevel/genelevel_counts_NTM_test_{minCount_expr}_{minSamples_expr}.tsv",
        TxIDs        = "resources/RNA/genelevel/genelevel_IDs_NTM_{minCount_expr}_{minSamples_expr}",
        tmm_fit      = "resources/RNA/genelevel/genelevel_tmm_fit_{minCount_expr}_{minSamples_expr}.rds"
    input:
        NT    = genelevel_raw("NT"),
        M     = genelevel_raw("M"),
        train = "resources/split/train_samples.txt",
        test  = "resources/split/test_samples.txt"
    conda:
        "../envs/r_get_transcriptome.yml"
    log:
        "logs/filter_counts_NTM/filter_counts_NTM_genelevel_{minCount_expr}_{minSamples_expr}.log"
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

rule batch_correct_rna_genelevel:   # batch_correct_rna on the gene counts (same script)
    output:
        counts_train = "resources/RNA/genelevel/genelevel_counts_NTM_train_combat_{minCount_expr}_{minSamples_expr}.tsv",
        counts_test  = "resources/RNA/genelevel/genelevel_counts_NTM_test_combat_{minCount_expr}_{minSamples_expr}.tsv",
        params_npz   = "resources/RNA/genelevel/genelevel_combat_params_{minCount_expr}_{minSamples_expr}.npz",
        plot_ds      = "results/plots/batch_correct/rna_genelevel_combat_by_dataset_{minCount_expr}_{minSamples_expr}.pdf",
        plot_cond    = "results/plots/batch_correct/rna_genelevel_combat_by_condition_{minCount_expr}_{minSamples_expr}.pdf"
    input:
        counts_train = "resources/RNA/genelevel/genelevel_counts_NTM_train_{minCount_expr}_{minSamples_expr}.tsv",
        counts_test  = "resources/RNA/genelevel/genelevel_counts_NTM_test_{minCount_expr}_{minSamples_expr}.tsv"
    params:
        sd_quantile_expr = config.get("sd_quantile_expr", 0)
    conda:
        "../envs/py_combat.yml"
    log:
        "logs/batch_correct/rna_genelevel_combat_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        PYTHONPATH=workflow/scripts/helpers python workflow/scripts/combat_correction_rna.py \
        --counts_train={input.counts_train} \
        --counts_test={input.counts_test} \
        --out_train={output.counts_train} \
        --out_test={output.counts_test} \
        --out_params={output.params_npz} \
        --out_plot_ds={output.plot_ds} \
        --out_plot_cond={output.plot_cond} \
        --sd_quantile={params.sd_quantile_expr} \
        --sd_row_prefix=ENST 2> {log}
        """
