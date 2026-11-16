STABLE_ONLY = int(config.get("stable_only", 1))   # 1: methyldriver / gimme use only the stably selected CpGs
if STABLE_ONLY:
    assert not CLUSTER, "stable_only = 1 needs cluster_cpgs: 0 (stability selection works on single CpGs)"

def sel_feats_input(wc):
    if STABLE_ONLY:
        return rules.stable_sel_feats.output[0]
    return "results/checks/model_selection" + TAG + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}"

rule selected_features_to_bed:
    output:
        bed_span   = "data/to_methyldriver" + TAG + "/selected_features_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.bed",
        bed_percpg = "data/to_methyldriver" + TAG + "/selected_features_percpg_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.bed"
    input:
        sel_feats = sel_feats_input
    params:
        sel_feats_dir = lambda wc, input: input.sel_feats if STABLE_ONLY else "results/sel_feats" + TAG,
        map_dir       = "data/clustering/clusters_maps" + TAG
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
