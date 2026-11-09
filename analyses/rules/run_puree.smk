rule run_puree:
    output:
        purity = "resources/purity/puree_purity_{sample_type}.tsv"
    input:
        matrix = "resources/RNA/puree_input_{sample_type}.tsv"
    conda:
        "../envs/py_puree.yml"
    log:
        "logs/puree/run_puree_{sample_type}.log"
    shell:
        """
        data_path=$(pwd)/{input.matrix}
        output_path=$(pwd)/{output.purity}
        log_path=$(pwd)/{log}
        cd resources/PUREE && python3 predict_purity.py \
        --data_path=$data_path \
        --output=$output_path \
        --gene_identifier_type=HGNC 2> $log_path
        """
