rule export_meth_metadata:
    output:
        coldata   = "resources/integrated_wgbs_hdf5/meth_colData.tsv",
        rowranges = "resources/integrated_wgbs_hdf5/meth_rowRanges.tsv"
    input:
        coldata   = "resources/integrated_wgbs_hdf5/integrated_wgbs_colData.rds",
        rowranges = "resources/integrated_wgbs_hdf5/integrated_wgbs_rowRanges.rds"
    conda:
        "../envs/dataset_integration.yml"
    log:
        "logs/batch_correct/export_meth_metadata.log"
    shell:
        """
        Rscript workflow/scripts/export_meth_metadata.R \
        --in_coldata={input.coldata} \
        --in_rowranges={input.rowranges} \
        --out_coldata={output.coldata} \
        --out_rowranges={output.rowranges} 2> {log}
        """