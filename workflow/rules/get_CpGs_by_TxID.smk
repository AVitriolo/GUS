rule get_CpGs_by_TxID:
    output:
        temp("data/{current_date}_{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}.TSS.distance_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_CpGs_by_{TxID}")
    input:
        "data/{current_date}_CpGs_coord_{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}.TSS.distance_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}.filtered.bed",
    shell:
        """
        grep {wildcards.TxID} {input} | cut -f 1 > {output} || echo "No matching CpGs found or {wildcards.TxID} is not in {input}"
        """