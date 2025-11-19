rule filter_coords_by_mappability:
    output:
        "resources/CpG_data/coord/{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.filtered.bed"
    input:
        mappability="resources/reference_data/mappability/{assembly_code}_k100.bismap.filtered.bedgraph",
        cpg_coord="resources/CpG_data/coord/{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.bed"
    conda:
        "/home/Prostate_GAS/workflow/envs/bedtools.yml"
    log:
        "logs/filter_coords_by_mappability/filter_coords_by_mappability_{assembly_code}_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.log"
    shell:
        """
	    bedtools intersect -u -a {input.cpg_coord} -b {input.mappability} > {output} 2> {log}
        """
