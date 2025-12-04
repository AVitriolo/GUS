rule plot_genome_tracks:
    output:
        "results/plots/tracks/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{TxID}_{importance}.pdf"
    input:
        tracks="results/genome_tracks/tracks/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{TxID}_{importance}.ini",
        TSS="resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.sorted.filtered.subset.bed"
    conda:
        "../envs/py_pyGenomeTracks.yml"
    log:
        "logs/plot_genome_tracks/plot_genome_tracks_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{TxID}_{importance}.log"
    shell:
        """
        pyGenomeTracks \
        --tracks {input.tracks} \
        --BED {input.TSS} \
        -out {output} 2> {log}
        """