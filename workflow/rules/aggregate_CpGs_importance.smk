rule aggregate_CpGs_importance:
    output:
        "results/genome_tracks/bed/aggr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{importance}.bed"
    input:
        check="results/checks/tracks/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{importance}"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/aggregate_CpGs_importance/aggregate_CpGs_importance_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{importance}.log"
    shell:
        """
        cat results/genome_tracks/bed/*_{wildcards.importance}.bed | bedtools sort -i - > {output} 2> {log}
        """