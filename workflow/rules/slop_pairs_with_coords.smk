rule slop_pairs_with_coords:
    output:
        "results/pairs/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}.slopped.bed"
    input:
        sorted_pairs_with_coords="results/pairs/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}.sorted.bed",
        chrom_sizes="resources/reference_data/{assembly_code}.chrom.clean.sizes"
    params:
        slop_distance_pairs_with_coords=config["slop_distances_pairs_with_coords"][0]
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/slop_pairs_with_coords/slop_pairs_with_coords_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}.log"
    shell:
        """
        bedtools slop -i {input.sorted_pairs_with_coords} -b {params.slop_distance_pairs_with_coords} -g {input.chrom_sizes} > {output} 2> {log}
        """