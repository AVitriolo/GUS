rule compute_pairwise_correlation:
    output:
        corr = temp("data/corr_NTM/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}")   # only merge_CpGs reads it (GBs per gene at +-1 Mb)
    input:
        rse   = "data/rse_NTM/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}",
        train = "resources/split/train_samples.txt"
    conda:
        "../envs/r_merge_CpGs.yml"
    priority: 20   # finish a gene (unpack -> corr -> clusters -> matrix -> model) before unpacking more: bounds temp RSEs on disk
    log:
        "logs/compute_pairwise_correlation/compute_pairwise_correlation_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{TxID}.log"
    shell:
        """
        Rscript workflow/scripts/compute_pairwise_correlation_split.R \
        --input_path_rse={input.rse} \
        --train_samples={input.train} \
        --output_path_corr={output.corr} \
        --tx_id={wildcards.TxID} 2> {log}
        """