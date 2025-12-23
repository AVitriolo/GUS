rule generate_genome_tracks:
    output:
        "results/genome_tracks/bed/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}_{importance}.bed"
    input:
        coords="results/genome_tracks/coords/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}",
        feat_imp="results/feat_imp/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}"
    conda:
        "../envs/r_CpGs.yml"
    log:
        "logs/generate_genome_tracks/generate_genome_tracks_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}_{importance}.log"
    shell:
        """
        Rscript workflow/scripts/generate_genome_tracks.R \
        --input_path_coords={input.coords} \
        --input_path_feat_imp={input.feat_imp} \
        --importance={wildcards.importance} \
        --output_path={output} 2> {log}
        """
        