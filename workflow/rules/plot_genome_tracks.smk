rule plot_genome_tracks:
    output:
        "results/plots/tracks/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{importance}"
    input:
        tracks="results/genome_tracks/tracks/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{importance}.ini",
        TSS="resources/reference_data/gencode/{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}.sorted.subset.slopped.bed"
    conda:
        "../envs/py_pyGenomeTracks.yml"
    log:
        "logs/plot_genome_tracks/plot_genome_tracks_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{comparison_type}_{minCount_expr}_{minSamples_expr}_{K}_{importance}.log"
    shell:
        """
        touch {output} && \
        pyGenomeTracks \
        --tracks {input.tracks} \
        --BED {input.TSS} \
        -out {output}.pdf 2> {log}
        """