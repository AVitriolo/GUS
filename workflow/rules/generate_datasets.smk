rule generate_datasets:
    output:
        "data/{current_date}_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_{comparison_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_xgbInput_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}"
    input:
        subsetted_beta="data/{current_date}_{assembly_code}_CpGs_beta_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_filtered_{comparison_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_subsetted_{TxID}",
        transcriptome="data/RNA/{current_date}_{sample_type}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{TxID}_transcriptome"
    shell:
        """
        Rscript workflow/scripts/generate_datasets.R \
        --input_transcriptome={input.transcriptome} \
        --input_subset_beta={input.subsetted_beta} \
        --TxID={wildcards.TxID} \
        --K={wildcards.K} \
        --output_path={output}
        """