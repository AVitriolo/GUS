rule selected_features_to_bed:
    output:
        bed_span   = "data/to_methyldriver/selected_features_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.bed",
        bed_percpg = "data/to_methyldriver/selected_features_percpg_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.bed"
    input:
        sel_check = "results/checks/model_selection/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}"
    params:
        sel_feats_dir = "results/sel_feats",
        map_dir       = "data/clustering/clusters_maps"
    conda:
        "../envs/r_merge_CpGs.yml"
    log:
        "logs/selected_features_to_bed/selected_features_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.log"
    shell:
        """
        Rscript workflow/scripts/selected_features_to_bed.R \
        --sel_feats_dir={params.sel_feats_dir} \
        --map_dir={params.map_dir} \
        --out_bed={output.bed_span} \
        --mode=span 2> {log}
        Rscript workflow/scripts/selected_features_to_bed.R \
        --sel_feats_dir={params.sel_feats_dir} \
        --map_dir={params.map_dir} \
        --out_bed={output.bed_percpg} \
        --mode=percpg 2>> {log}
        """