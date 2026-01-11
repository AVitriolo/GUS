rule get_TxID_specific_gtf:
    output:
        "resources/reference_data/gencode/gencode.v{gencode_version}_{TxID}.annotation.gtf.gz"
    input:
        "resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gtf.gz",
    log:
        "logs/get_TxID_specific_gtf/get_TxID_specific_gtf_v{gencode_version}_{TxID}.log"
    shell:
        """
        zcat {input} | grep -Fw {wildcards.TxID} | gzip > {output}
        """