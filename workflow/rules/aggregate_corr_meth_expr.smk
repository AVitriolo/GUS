rule aggregate_corr_meth_expr:
    output:
        "results/corr_meth_expr/aggr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}"
    input:
        checks="results/checks/corr_meth_expr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}"
    log:
        "logs/corr_meth_expr/aggregate_corr_meth_expr_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}.log"
    shell:
        """
        for f in results/corr_meth_expr/* ; do
            if [ -s "$f" ] && [ -f "$f" ]; then
                enst=$(basename "$f" | awk -F'_' '{{print $NF}}')
                awk -v id="$enst" 'FNR>1 && NF {{print $0 "\t" id}}' "$f" >> {output}

            fi
        done
        """
