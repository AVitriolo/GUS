rule methyldriver:                      # runs PER comparison, cheap
    output:
        results    = "results/methyldriver/{comparison}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.rds",
        volcano    = "results/methyldriver/volcano_{comparison}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.pdf",
        full_table = "results/methyldriver/full_{comparison}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.tsv"
    input:
        neighbors = "results/methyldriver/neighbors_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.rds",
        matrix    = "results/methyldriver/featmatrix_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.rds",
        meth_diff = "data/to_methyldriver/{comparison}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.tsv"
    params:
        conversion_table = config["conversion_table"], q = 0.05
    conda: "../envs/r_methyldriver.yml"
    log: "logs/methyldriver/{comparison}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.log"
    shell:
        """
        Rscript workflow/scripts/run_methyldriver.R \
        --input_neighbors={input.neighbors} --input_matrix={input.matrix} \
        --input_meth_diff={input.meth_diff} --conversion_table={params.conversion_table} \
        --output_results={output.results} --output_volcano={output.volcano} \
        --output_full_table={output.full_table} --q={params.q} > {log} 2>&1
        """