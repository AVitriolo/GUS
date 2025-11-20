rule garbage_collector:
    output:
        "results/checks/garbage_collector/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{importance}"
    input:
        hyper="results/hyper/aggr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}.txt",
        performance="results/performance/aggr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}.txt",
        tracks="results/genome_tracks/bed/aggr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{importance}.bed"
    log:
        "logs/garbage_collector/garbage_collector_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{importance}.log"
    shell:
        """
        pattern_hyper=results/hyper/{wildcards.assembly_code}_{wildcards.sample_type}_{wildcards.leftCount_beta}_{wildcards.rightCount_beta}_{wildcards.minSamples_beta}_{wildcards.minCov}_v{wildcards.gencode_version}_{wildcards.tss_subset}_{wildcards.distance}_{wildcards.min_CpG}_{wildcards.comparison_type}_{wildcards.minCount_expr}_{wildcards.minSamples_expr}_{wildcards.K}_*;
        pattern_performance=results/performance/{wildcards.assembly_code}_{wildcards.sample_type}_{wildcards.leftCount_beta}_{wildcards.rightCount_beta}_{wildcards.minSamples_beta}_{wildcards.minCov}_v{wildcards.gencode_version}_{wildcards.tss_subset}_{wildcards.distance}_{wildcards.min_CpG}_{wildcards.comparison_type}_{wildcards.minCount_expr}_{wildcards.minSamples_expr}_{wildcards.K}_*
        pattern_bed=results/genome_tracks/bed/{wildcards.assembly_code}_{wildcards.sample_type}_{wildcards.leftCount_beta}_{wildcards.rightCount_beta}_{wildcards.minSamples_beta}_{wildcards.minCov}_v{wildcards.gencode_version}_{wildcards.tss_subset}_{wildcards.distance}_{wildcards.min_CpG}_{wildcards.comparison_type}_{wildcards.minCount_expr}_{wildcards.minSamples_expr}_{wildcards.K}_*_{wildcards.importance}.bed
        touch {output}
        for file in $pattern_hyper; do rm $file; done
        for file in $pattern_performance; do rm $file; done
        for file in $pattern_bed; do rm $file; done
        """