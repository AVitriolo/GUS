rule get_GRN:
    output:
        "results/GRN/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}"
    input:
        path_gimme="results/gimme/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}.processed",
        path_converter=config["converter_location"][0],
        path_slopped="results/pairs/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}.slopped.bed",
        path_gff="resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gff3.gz",
        path_expressed_TxIDs="resources/TxIDs/TxIDs_{minCount_expr}_{minSamples_expr}"
    conda:
        "../envs/r_GRN.yml"
    log:
        "logs/get_GRN/get_GRN_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}.log"
    shell:
        """
        Rscript workflow/scripts/get_GRN.R \
        --input_path_gimme={input.path_gimme} \
        --input_path_converter={input.path_converter} \
        --input_path_slopped={input.path_slopped} \
        --input_path_gff={input.path_gff} \
        --input_path_expressed_TxIDs={input.path_expressed_TxIDs} \
        --output_path={output} 2> {log}
        """