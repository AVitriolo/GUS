def puree_counts_input(wildcards):
    if wildcards.sample_type == "M":
        return config["RNA_counts_locations_M"][0]
    return config["RNA_counts_locations_NT"][0]

def puree_sample_prefix_arg(wildcards):
    if wildcards.sample_type == "M":
        return ""
    return f"--sample_prefix={wildcards.sample_type}"

rule prepare_puree_input:
    output:
        matrix = "resources/RNA/puree_input_{sample_type}.tsv"
    input:
        counts = puree_counts_input,
        conversion_table = config["conversion_table"]
    params:
        sample_prefix_arg = puree_sample_prefix_arg
    conda:
        "../envs/py_puree.yml"
    log:
        "logs/puree/prepare_puree_input_{sample_type}.log"
    shell:
        """
        python3 analyses/scripts/prepare_puree_input.py \
        --counts={input.counts} \
        --conversion_table={input.conversion_table} \
        {params.sample_prefix_arg} \
        --output={output.matrix} 2> {log}
        """
