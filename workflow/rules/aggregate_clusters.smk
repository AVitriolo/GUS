rule aggregate_clusters:
    output:
        "results/corr_meth_expr/bed/clustered/aggr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}"
    input:
        checks="results/checks/cluster/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}"
    log:
        "logs/cluster_corr_meth_expr/aggr/aggr_clusters_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}.log"
    shell:
        """
        for f in results/corr_meth_expr/bed/clustered/* ; do
            [ -f "$f" ] || continue
            [ -s "$f" ] || continue
            awk 'NR>1 {{print $1,$2,$3,$4,$7}}' "$f" >> {output}
        done
        """