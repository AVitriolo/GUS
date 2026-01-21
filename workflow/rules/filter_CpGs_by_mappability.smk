rule filter_CpGs_by_mappability:
    output:
        CpGs="resources/CpG_data/coord/{assembly_code}.bed"
    input:
        dir=config["CpGs_locations"][0],
        mappability="resources/reference_data/mappability/{assembly_code}_k100.bismap.filtered.bedgraph"
    conda:
        "../envs/r_CpGs.yml"
    log:
        "logs/filter_CpGs_by_mappability/filter_CpGs_by_mappability_{assembly_code}.log"
    shell:
        """
         Rscript workflow/scripts/filter_CpGs_by_mappability.R \
        --input_dir={input.dir} \
        --input_path_mappability={input.mappability} \
        --output_path={output}
        """