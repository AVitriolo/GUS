rule filter_mappability:
    output:
        "resources/reference_data/mappability/{assembly_code}_k100.bismap.filtered.bedgraph"
    input:
        "resources/reference_data/mappability/{assembly_code}_k100.bismap.bedgraph"
    conda:
        "../envs/r_filter_CpGs_TSSs.yml"
    log:
        "logs/filter_mappability/filter_mappability_{assembly_code}.log"
    shell:
        """
        Rscript workflow/scripts/filter_mappability.R \
        --input_path={input} \
        --output_path={output}
	    """