rule get_TxIDs_to_exclude:
    output:
        "resources/TxIDs/TxIDs_to_exclude"
    log:
        "logs/get_TxIDs_to_exclude/get_TxIDs_to_exclude.log"
    shell:
        """
        touch {output} 2> {log}
        """