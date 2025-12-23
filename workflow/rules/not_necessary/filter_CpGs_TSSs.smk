rule filter_CpGs_TSSs:
    output:
        distance_filtered="data/coord/{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.distance.filtered.bed",
        plot_1="results/plots/filtering/nCpG_byTx/{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.pdf",
        plot_2="results/plots/filtering/nCpG_byTx_distance/{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.pdf"
    input:
        "data/coord/{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.distance.bed"
    conda:
        "../envs/r_filter_CpGs_TSSs.yml"
    log:
        "logs/filter_CpGs_TSSs/filter_CpGs_TSSs_{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.log"
    shell:
        """
        Rscript workflow/scripts/filter_CpGs_TSSs.R \
        --input_path={input} \
        --distance_cutoff={wildcards.distance} \
        --min_CpG={wildcards.min_CpG} \
        --output_distance_filtered_path={output.distance_filtered} \
        --output_plot_1_path={output.plot_1} \
        --output_plot_2_path={output.plot_2} 2> {log}
        """