rule cluster_delta_methylation:
    output:
        tsv = "data/to_methyldriver/{comparison}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.tsv"
    input:
        h5  = integrated_h5,
        bed = "data/to_methyldriver/selected_features_percpg_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.bed"
    params:
        cond_1 = lambda wc: config["delta_comparisons"][wc.comparison]["cond_1"],
        cond_2 = lambda wc: config["delta_comparisons"][wc.comparison]["cond_2"],
        site_2 = lambda wc: (
            f"--site_2={config['delta_comparisons'][wc.comparison]['site_2']}"
            if "site_2" in config["delta_comparisons"][wc.comparison] else ""
        )
    conda:
        "../envs/r_merge_CpGs.yml"
    log:
        "logs/to_methyldriver/{comparison}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.log"
    shell:
        """
        Rscript workflow/scripts/cluster_delta_methylation.R \
        --input_h5={input.h5} \
        --input_bed={input.bed} \
        --cond_1={params.cond_1} \
        --cond_2={params.cond_2} {params.site_2} \
        --output={output.tsv} 2> {log}
        """