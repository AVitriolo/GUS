rule clean_chrom_sizes:
    output:
        "resources/reference_data/{assembly_code}.chrom.clean.sizes"
    input:
        "resources/reference_data/{assembly_code}.chrom.sizes"
    shell:
        """
	    cat {input} | grep -v "_" > {output} 
        """