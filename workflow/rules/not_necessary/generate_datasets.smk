rule generate_datasets:
    output:
        xgb = "data/xgb/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}",
        corr_table = "data/corr/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}"
    input:
        subsetted_beta="data/beta/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_{comparison_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{TxID}",
        transcriptome="data/RNA/{sample_type}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{TxID}"
    conda:
        "../envs/r_CpGs.yml"
    log:
        "logs/generate_datasets/generate_datasets_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/generate_datasets.R \
        --input_transcriptome={input.transcriptome} \
        --input_subset_beta={input.subsetted_beta} \
        --TxID={wildcards.TxID} \
        --K={wildcards.K} \
        --output_path_xgb={output.xgb} \
        --output_path_corr_table={output.corr_table} 2> {log}
        """