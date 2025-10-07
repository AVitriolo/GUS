rule get_chrom_sizes:
    output:
        "ReferenceData/{assembly_code}.chrom.clean.sizes"
    input:
        "ReferenceData/{assembly_code}.chrom.sizes"
    shell:
        """
	    cat {input} | grep -v "_" > {output} 
        """