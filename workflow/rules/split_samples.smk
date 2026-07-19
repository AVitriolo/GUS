rule split_samples:
    output:
        train   = "resources/split/train_samples.txt",
        test    = "resources/split/test_samples.txt",
        kept    = "resources/split/kept_samples.txt",
        paired  = "resources/split/paired_samples.txt",
        summary = "resources/split/split_summary.tsv",
        plot    = "results/plots/split/split_composition.pdf"
    input:
        NT   = config["RNA_counts_locations_NT"][0],
        M    = config["RNA_counts_locations_M"][0],
        wgbs = "resources/split/wgbs_samples.txt"
    params:
        test_size   = config["test_sizes"][0],
        seed        = config["split_seeds"][0],
        min_libsize = config["min_libsizes"][0],
        stratify    = "--stratify" if config["stratify_splits"][0] else "",
        outdir      = "resources/split"
    conda:
        "../envs/py_split.yml"
    log:
        "logs/split_samples/split_samples.log"
    shell:
        """
        python workflow/scripts/train_test_split.py \
        --input_expr_NT={input.NT} \
        --input_expr_M={input.M} \
        --input_wgbs_samples={input.wgbs} \
        --min_libsize={params.min_libsize} \
        --test_size={params.test_size} \
        --seed={params.seed} \
        --outdir={params.outdir} \
        --out_plot={output.plot} \
        {params.stratify} 2> {log}
        """
