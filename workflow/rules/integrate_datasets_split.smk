rule integrate_datasets:
    output:
        h5        = "resources/integrated_wgbs_hdf5/integrated_wgbs.h5",
        rowranges = "resources/integrated_wgbs_hdf5/integrated_wgbs_rowRanges.rds",
        coldata   = "resources/integrated_wgbs_hdf5/integrated_wgbs_colData.rds",
        fit       = "resources/integrated_wgbs_hdf5/integrated_wgbs_train_fit.rds"
    input:
        cpgea = config["cpgea_wgbs_location"][0],
        mcrpc = config["mcrpc_wgbs_location"][0],
        train = "resources/split/train_samples.txt",
        test  = "resources/split/test_samples.txt"
    conda:
        "../envs/dataset_integration.yml"
    log:
        "logs/integrate_datasets/integrate_datasets.log"
    shell:
        """
        Rscript workflow/scripts/dataset_integration_split.R \
        --cpgea={input.cpgea} \
        --mcrpc={input.mcrpc} \
        --train_samples={input.train} \
        --test_samples={input.test} \
        --out_h5={output.h5} \
        --out_rowranges={output.rowranges} \
        --out_coldata={output.coldata} \
        --out_fit={output.fit} 2> {log}
        """