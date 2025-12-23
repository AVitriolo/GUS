rule plot_genome_tracks:
    output:
        "results/plots/tracks/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}_{importance}.pdf"
    input:
        tracks="results/genome_tracks/tracks/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}_{importance}.ini",
        TSS="resources/reference_data/gencode/slopped/{assembly_code}_v{gencode_version}_{tss_subset}_{minCount_expr}_{minSamples_expr}_{distance}_{TxID}.slopped.bed"
    conda:
        "../envs/py_pyGenomeTracks.yml"
    log:
        "logs/plot_genome_tracks/plot_genome_tracks_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}_{importance}.log"
    shell:
        """
        temp_file="results/plots/tracks/temp_{wildcards.TxID}.pdf"
        temp_pattern="results/plots/tracks/temp_{wildcards.TxID}_*.pdf"
        pyGenomeTracks --tracks {input.tracks} --BED {input.TSS} -out $temp_file 2> {log}
        mv $temp_pattern {output}
        """