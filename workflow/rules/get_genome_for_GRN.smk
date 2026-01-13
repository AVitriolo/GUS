rule get_genome_for_GRN:
    output:
        "results/checks/GRN/{assembly_code}"
    conda:
        "../envs/gimme.yml"
    log:
        "logs/get_genome_for_GRN/get_genome_for_GRN_{assembly_code}.log"
    shell:
        """
        touch {output}
        genomepy install GRCh38.p14 --provider Ensembl --mask soft 2>{log}
        """