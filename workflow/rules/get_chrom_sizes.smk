rule get_chrom_sizes:
    output:
	    "resources/reference_data/{assembly_code}.chrom.sizes"
    log:
        "logs/get_chrom_sizes/get_chrom_sizes_{assembly_code}.log"
    shell:
        """
	    wget -qP resources/reference_data/ https://hgdownload.cse.ucsc.edu/goldenpath/{wildcards.assembly_code}/bigZips/{wildcards.assembly_code}.chrom.sizes 2> {log}
        """