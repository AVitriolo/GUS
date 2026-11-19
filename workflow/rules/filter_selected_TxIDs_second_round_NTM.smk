checkpoint filter_selected_TxIDs_second_round_NTM:
    # streaming split of the CpG x window intersect into one candidate-CpG file per TxID (split_topK_by_TxID.sh);
    # the old R script loaded the whole table (33 GB for +-1 Mb) and was killed for memory
    output:
        ff_selected_TxIDs  = "resources/TxIDs/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_ff_selected_TxIDs",
        split_dir          = directory("resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_split"),
        cpgs_per_txid      = "results/plots/filtering/nCpGs/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.tsv",
        distance_hist      = "results/plots/filtering/distance/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.tsv",
        path_num_CpGs_plot = "results/plots/filtering/nCpGs/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.pdf",
        path_dist_filt_plot= "results/plots/filtering/distance/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.pdf"
    input:
        path_CpGs = "resources/CpG_data/top_K/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.bed"
    conda:
        "../envs/r_xgb_input.yml"
    log:
        "logs/filter_selected_TxIDs_second_round/filter_selected_TxIDs_second_round_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}.log"
    shell:
        """
        bash workflow/scripts/split_topK_by_TxID.sh {input.path_CpGs} {output.split_dir} \
            {output.ff_selected_TxIDs} {output.cpgs_per_txid} {output.distance_hist} > {log} 2>&1
        Rscript workflow/scripts/plot_topK_filter_stats.R {output.cpgs_per_txid} {output.distance_hist} \
            {output.path_num_CpGs_plot} {output.path_dist_filt_plot} >> {log} 2>&1
        """
