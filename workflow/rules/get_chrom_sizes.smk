rule get_chrom_sizes:
    output:
	    "ReferenceData/{assembly_code}.chrom.sizes"
    shell:
        """
	    wget -qP ReferenceData/ https://hgdownload.cse.ucsc.edu/goldenpath/{wildcards.assembly_code}/bigZips/{wildcards.assembly_code}.chrom.sizes 
        """