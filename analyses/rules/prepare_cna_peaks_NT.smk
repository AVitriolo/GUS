rule prepare_cna_peaks_NT:
    output:
        peaks = "resources/CNAs_NT/cna_peaks_NT.tsv"
    input:
        xlsx = config["cna_peaks_xlsx"]
    params:
        sheet = config["cna_peaks_sheet"]
    conda:
        "../envs/r_camdac.yml"
    log:
        "logs/cna_correction/prepare_cna_peaks_NT.log"
    shell:
        """
        Rscript analyses/scripts/prepare_cna_peaks_NT.R \
        --input_xlsx={input.xlsx} \
        --sheet="{params.sheet}" \
        --output_path={output.peaks} 2> {log}
        """
