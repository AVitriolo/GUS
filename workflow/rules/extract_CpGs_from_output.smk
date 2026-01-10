rule extract_CpGs_from_output:
    output:
        coords="results/genome_tracks/coords/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}"
    input:
        coords="resources/composite/{assembly_code}_{sample_type}.bed",
        patterns="results/sel_feats/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}"
    log:
        "logs/extract_CpGs_from_output/extract_CpGs_from_output_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}.log"
    shell:
        """
        grep -Fw -f {input.patterns} {input.coords} > {output.coords} 2> {log} || echo "No selected CpGs found for {wildcards.TxID} in {input.coords} >> {log}"
        """
