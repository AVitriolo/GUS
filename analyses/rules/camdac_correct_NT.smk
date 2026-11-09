rule camdac_correct_NT:
    output:
        h5        = "resources/integrated_wgbs_hdf5/integrated_wgbs_cna.h5",
        coldata   = "resources/integrated_wgbs_hdf5/integrated_wgbs_cna_colData.rds",
        rowranges = "resources/integrated_wgbs_hdf5/integrated_wgbs_cna_rowRanges.rds",
        summary   = "results/cna_correction/camdac_NT_summary.tsv"
    input:
        h5        = "resources/integrated_wgbs_hdf5/integrated_wgbs.h5",
        coldata   = "resources/integrated_wgbs_hdf5/integrated_wgbs_colData.rds",
        rowranges = "resources/integrated_wgbs_hdf5/integrated_wgbs_rowRanges.rds",
        peaks     = "resources/CNAs_NT/cna_peaks_NT.tsv",
        purity    = config["cna_purity_T"]
    params:
        min_cn = config["cna_min_cn"]
    conda:
        "../envs/r_camdac.yml"
    log:
        "logs/cna_correction/camdac_correct_NT.log"
    shell:
        """
        # the correction rewrites only the rows inside peaks, in place on a copy
        cp {input.h5} {output.h5}
        cp {input.coldata}   {output.coldata}
        cp {input.rowranges} {output.rowranges}
        Rscript analyses/scripts/camdac_correct_NT.R \
        --h5={output.h5} \
        --rowranges={input.rowranges} \
        --coldata={input.coldata} \
        --peaks={input.peaks} \
        --purity={input.purity} \
        --min_cn={params.min_cn} \
        --output_summary={output.summary} 2> {log}
        """
