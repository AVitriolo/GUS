rule slop_TSSs:
    output:
        "resources/reference_data/gencode/{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.sorted.subset.slopped.bed"
    input:
        TSS="resources/reference_data/gencode/{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.sorted.subset.bed",
        genome="resources/reference_data/{assembly_code}.chrom.clean.sizes"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/slop_TSSs/slop_TSSs_{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.log"
    shell:
        """
        bedtools slop -i {input.TSS} -g {input.genome} -b {wildcards.distance} > {output} 2> {log}
        """