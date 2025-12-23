rule generate_ini_file:
    output:
        "results/genome_tracks/tracks/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}_{importance}.ini"
    input:
        gtf="resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gtf.gz",
        CpG_importance="results/genome_tracks/bed/{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}_{importance}.bed"
    log:
        "logs/generate_ini_file/generate_ini_file_{assembly_code}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_{minCov}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_{minCount_expr}_{minSamples_expr}_{K_closest}_{test_size}_{CV}_{TxID}_{importance}.log"
    shell:
        r"""
cat << EOF > {output} 2> {log}
[x-axis]
where = top

[cpg_importance]
file = {input.CpG_importance}
height = 6
color = #1f77b4
file_type = bedgraph
min_value = 0
max_value = auto

[spacer]
height = 0.5

[genes]
file = {input.gtf}
height = 6
merge_transcripts = false
style = UCSC
prefered_name = transcript_id
fontsize = 12
file_type = gtf
arrow_length = 20
EOF
        """