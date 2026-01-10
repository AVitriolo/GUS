rule get_composite_coordinates:
    output:
        path_composite="resources/composite/{assembly_code}_{sample_type}.bed"
    input:
        TMR_file="resources/TMRs/{sample_type}_prostate_tmrs.bed",
        CpG_sorted_coord="resources/CpG_data/coord/{assembly_code}.sorted.bed"
    log:
        "logs/get_composite_coordinates/get_composite_coordinates_{assembly_code}_{sample_type}.log"
    shell:
        """
        cat <(cut -f 1-4 {input.CpG_sorted_coord}) <(cut -f 1-4 {input.TMR_file}) > {output.path_composite} 2> {log}
        """