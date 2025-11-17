rule clean_chrom_sizes:
    output:
        "resources/reference_data/{assembly_code}.chrom.clean.sizes"
    input:
        "resources/reference_data/{assembly_code}.chrom.sizes"
    log:
        "logs/clean_chrom_sizes/clean_chrom_sizes_{assembly_code}.log"
    shell:
        """
	    cat {input} | grep -v "_" > {output} 2> {log}
        """