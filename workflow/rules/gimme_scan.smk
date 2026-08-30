rule gimme_scan:
    output:
        scan = "results/gimme/scan_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.tsv"
    input:
        bed         = "data/to_methyldriver/selected_features_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.bed",
        chrom_sizes = "resources/reference_data/hg38.chrom.sizes"
    params:
        genome  = "GRCh38.p14",
        min_len = 20,
        pfm     = "",
        subset  = lambda wc: config.get("gimme_subset", "")   # "" = all regions
    conda:
        "../envs/gimme.yml"
    log:
        "logs/gimme/scan_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.log"
    shell:
        """
        bash workflow/scripts/gimme_scan.sh \
        {input.bed} {params.genome} {output.scan} {input.chrom_sizes} \
        {params.min_len} "{params.pfm}" "{params.subset}" 2> {log}
        """
