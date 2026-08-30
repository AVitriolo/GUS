rule methyldriver:
    output:
        results    = "results/methyldriver/{comparison}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.rds",
        volcano    = "results/methyldriver/volcano_{comparison}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.pdf",
        full_table = "results/methyldriver/full_{comparison}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.tsv"
    input:
        bed        = "data/to_methyldriver/selected_features_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.bed",
        meth_diff  = "data/to_methyldriver/{comparison}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.tsv"
    params:
        bigwig_dir       = config["encode_bigwig_dir"],
        conversion_table = config["conversion_table"],
        exclude_files    = "ENCFF582VMO.bigWig,ENCFF582VMO.bigWig.1",
        clusters_to_keep = lambda wc: (
            f"--clusters_to_keep={config['clusters_to_keep']}"
            if config.get("clusters_to_keep") else ""
        ),
        n      = 100,
        n_pcs  = 20,
        q      = 0.05
    threads: 8
    conda:
        "../envs/r_methyldriver.yml"
    log:
        "logs/methyldriver/{comparison}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.log"
    shell:
        """
        Rscript workflow/scripts/run_methyldriver.R \
        --input_bed={input.bed} \
        --input_meth_diff={input.meth_diff} \
        --bigwig_dir={params.bigwig_dir} \
        --conversion_table={params.conversion_table} \
        --exclude_files={params.exclude_files} \
        {params.clusters_to_keep} \
        --output_results={output.results} \
        --output_volcano={output.volcano} \
        --output_full_table={output.full_table} \
        --num_cores={threads} \
        --n={params.n} --n_pcs={params.n_pcs} --q={params.q} 2> {log}
        """
