rule sort_CpGs:
    output:
        "resources/reference_data/"
    input:
        "resources/reference_data/"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
	     sortBed -i {input} > {output}
        """