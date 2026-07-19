rule export_wgbs_samples:
    output:
        samples = "resources/split/wgbs_samples.txt"
    input:
        cpgea = config["cpgea_wgbs_location"][0],
        mcrpc = config["mcrpc_wgbs_location"][0]
    conda:
        "../envs/dataset_integration.yml"
    log:
        "logs/split_samples/export_wgbs_samples.log"
    shell:
        """
        Rscript workflow/scripts/export_wgbs_samples.R \
        --cpgea={input.cpgea} \
        --mcrpc={input.mcrpc} \
        --output_samples={output.samples} 2> {log}
        """