rule get_gencode_annotation:
    output:
	    "resources/reference_data/gencode.v{gencode_version}.annotation.gff3.gz"
    shell:
        """
        wget -qP resources/reference_data/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_version}/gencode.v{wildcards.gencode_version}.annotation.gff3.gz
        """