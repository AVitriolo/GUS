rule sort_CpGs:
    output:
        "resources/CpG_data/coord/{assembly_code}.sorted.bed"
    input:
        "resources/CpG_data/coord/{assembly_code}.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/sort_CpGs/sort_CpGs_{assembly_code}.log"
    shell:
        """
	     bedtools sort -i {input} > {output} 2> {log}
        """
        