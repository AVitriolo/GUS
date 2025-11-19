rule get_mappability:
    output:
        "resources/reference_data/mappability/{assembly_code}_k100.bismap.bedgraph"
    log:
        "logs/get_mappability/get_mappability_{assembly_code}.log"
    shell:
        """
        wget -qO- https://bismap.hoffmanlab.org/raw/{wildcards.assembly_code}/k100.bismap.bedgraph.gz | gunzip > {output} 2>> {log}
	    """
