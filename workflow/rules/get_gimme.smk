rule gimme:
    output:
        "results/gimme/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}"
    input:
        pairs="results/pairs/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}.prepared.bed",
        check="results/checks/GRN/{assembly_code}"
    conda:
        "../envs/gimme.yml"
    log:
        "logs/gimme/gimme_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}.log"
    shell:
        """
        gimme scan {input.pairs} -g GRCh38.p14 > {output} 2> {log}
        """