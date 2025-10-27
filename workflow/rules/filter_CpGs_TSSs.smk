rule filter_CpGs_TSSs:
    output:
        distance_filtered="data/{current_date}_CpGs_coord_{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}.TSS.distance_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}.filtered.bed",
        plot_1="results/plots/filtering/{current_date}_{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_nCpG_byTx_boxLog10.pdf",
        plot_2="results/plots/filtering/{current_date}_{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_nCpG_byTx_distance.pdf"
    input:
        "data/{current_date}_CpGs_coord_{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}.TSS.distance.bed"
    conda:
        "../envs/r_filter_CpGs_TSSs.yml"
    shell:
        """
        Rscript workflow/scripts/filter_CpGs_TSSs.R \
        --input_path={input} \
        --distance_cutoff={wildcards.distance} \
        --min_CpG={wildcards.min_CpG} \
        --output_distance_filtered_path={output.distance_filtered} \
        --output_plot_1_path={output.plot_1} \
        --output_plot_2_path={output.plot_2}
        """