rule subset_beta_matrix:
    output:
        temp("data/{current_date}_{assembly_code}_CpGs_beta_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_filtered_{comparison_type}_v{gencode_version}_{tss_subset}_{distance}_{min_CpG}_subsetted_{TxID}")
    input:
        patterns="data/{current_date}_{comparison_type}_{assembly_code}_v{gencode_version}_{tss_subset}.TSS.distance_{distance}_{min_CpG}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}_CpGs_by_{TxID}",
        betas="resources/CpG_data/{current_date}_{assembly_code}_CpGs_beta_{comparison_type}_{sample_type}_{leftCount_beta}_{rightCount_beta}_{minSamples_beta}.txt"
    shell:
        """
        head -n 1 {input.betas} > {output}
        grep -F -f {input.patterns} {input.betas} >> {output}
        """
