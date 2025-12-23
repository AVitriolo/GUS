rule get_K_closest_CpGs:
    output:
        "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}_{K_closest}.bed"
    input:
        CpGs="resources/CpG_data/coord/{assembly_code}.sorted.bed",
        filtered_TSSs="resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}.sorted.filtered.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/get_K_closest_CpGs/get_K_closest_CpGs_{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}_{K_closest}.log"
    shell:
        """
        bedtools closest -a {input.filtered_TSSs} -b {input.CpGs} -t all -k {wildcards.K_closest} -d | awk -F '\t' 'BEGIN{{OFS="\t"}} {{print $5, $6, $7, $8, $4, $9}}' > {output} 2>{log}
        """



