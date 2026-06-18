rule intersect_CpGs_TSSs:
    output:
        "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.bed"
    input:
        CpGs   = "resources/CpG_data/coord/{assembly_code}.sorted.bed",
        window = "resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}_{distance}.window.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/intersect_CpGs_TSSs/intersect_CpGs_TSSs_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        bedtools intersect -a {input.CpGs} -b {input.window} -wa -wb | \
        awk -F'\t' 'BEGIN{{OFS="\t"}} {{
            tss = $6 + {wildcards.distance};
            dist = $2 - tss;
            if (dist < 0) dist = -dist;
            print $1, $2, $3, $4, $8, dist
        }}' > {output} 2> {log}
        """
