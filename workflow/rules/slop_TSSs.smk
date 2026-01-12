rule slop_TSSs:
    output:
        "resources/reference_data/gencode/slopped/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}_{distance}_{TxID}.slopped.bed"
    input:
        TSS="resources/reference_data/gencode/TSS_to_plot/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}_{TxID}.bed",
        genome="resources/reference_data/{assembly_code}.chrom.clean.sizes"
    params:
        slop_distance_TSS=config["slop_distances_TSS"][0]
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/slop_TSSs/slop_TSSs_{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}_{distance}_{TxID}.log"
    shell:
        """
        bedtools slop -i {input.TSS} -g {input.genome} -b {params.slop_distance} > {output} 2> {log}
        """