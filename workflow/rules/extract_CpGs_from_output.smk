rule extract_CpGs_from_output:
    output:
        coords="results/genome_tracks/coords/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}"
    input:
        coords="resources/CpG_data/coord/{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.sorted.bed",
        patterns="results/sel_feats/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}",
        check="results/checks/train/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}"
    log:
        "logs/extract_CpGs_from_output/extract_CpGs_from_output_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}.log"
    shell:
        """
        grep -F -f {input.patterns} {input.coords} > {output.coords} 2> {log} || echo "No selected CpGs found for {wildcards.TxID} in {input.coords} >> {log}"
        """
