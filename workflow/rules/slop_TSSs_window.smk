rule slop_TSSs_window:
    output:
        "resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}_{distance}.window.bed"
    input:
        TSSs       = "resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}.sorted.filtered.bed",
        chrom_sizes = "resources/reference_data/{assembly_code}.chrom.clean.sizes"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/slop_TSSs_window/slop_TSSs_window_{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}_{distance}.log"
    shell:
        """
        bedtools slop -i {input.TSSs} -g {input.chrom_sizes} -b {wildcards.distance} > {output} 2> {log}
        """
