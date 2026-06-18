rule integrate_wgbs:
    input:
        cpgea = ancient("resources/cpgea_wgbs_with_coverage_hg38"),
        mcrpc = ancient("resources/mcrpc_wgbs_hg38")
    output:
        h5        = "resources/integrated_wgbs_hdf5/integrated_wgbs.h5",
        rowranges = "resources/integrated_wgbs_hdf5/integrated_wgbs_rowRanges.rds",
        coldata   = "resources/integrated_wgbs_hdf5/integrated_wgbs_colData.rds"
    conda:
        "../envs/dataset_integration.yml"
    log:
        "logs/integrate_wgbs/integrate_wgbs.log"
    shell:
        """
        Rscript workflow/scripts/dataset_integration.R \
        --cpgea {input.cpgea} \
        --mcrpc {input.mcrpc} \
        --out_h5 {output.h5} \
        --out_rowranges {output.rowranges} \
        --out_coldata {output.coldata} 2> {log}
        """
