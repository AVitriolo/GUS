rule process_gimme:
    output:
        "results/gimme/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}.processed"
    input:
        "results/gimme/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}"
    log:
        "logs/process_gimme/process_gimme_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}.log"
    shell:
        """
        paste <(tail -n+6 {input} | cut -d " " -f 3 | sed 's/"//g') \
              <(cut -f 1,6 <(tail -n+6 {input} | cut -d " " -f 2)) > {output} 2>{log}
        """