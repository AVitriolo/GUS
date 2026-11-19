# ---- stability selection (Shah & Samworth complementary pairs, PFER <= 1) on the models' selected CpGs ----
# Same chain as the genestd run (analyses/scripts, run by hand before):
#   1. stability_selection_kbest.py   (sharded) the models' SelectKBest(f_regression, gene's k) on 2 x n_pairs half-samples
#   2. stability_selection_summary.R  per-gene stability, per-CpG frequency (cpg_freq.tsv)
#   3. stability_regions_corr.py      (sharded) co-methylated regions (adjacent r >= corr, <= max_gap bp)
#   4. stability_selection_regions.R  region-level frequency + PFER cutoff (region_freq_corr<r>.tsv)
#   5. stable_sel_feats.py            each model's selected CpGs lying in a stable region -> sel_feats format
# selected_features_to_bed (-> methyldriver, gimme) reads 5. when stable_only = 1.
# cluster_cpgs: 0 -> the chain above, on single CpGs (as the genestd run).
# cluster_cpgs: 1 -> candidates are the model's cluster features: step 1 = stability_selection_features.py on the
#                    model's own training matrix, steps 3-4 are skipped (a cluster is already a co-methylated
#                    region), step 5 = stable_sel_feats_features.py (features with frequency >= the PFER cutoff).

STAB_DIR    = "results/stability_selection" + TAG + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}"
STAB_SHARDS = int(config.get("stability_shards", 24))
STAB_PAIRS  = int(config.get("stability_n_pairs", 100))
STAB_PFER   = float(config.get("stability_pfer", 1))
STAB_CORR   = str(config.get("stability_corr", 0.5))
STAB_GAP    = int(config.get("stability_max_gap", 300))

def stab_prefix(wc):   # prefix of the per-gene candidate BEDs (resources/CpG_data/top_K/<prefix>_<TxID>.bed)
    return f"{wc.assembly_code}_v{wc.gencode_version}_{wc.tss_subset}_{wc.distance}_{wc.minCount_expr}_{wc.minSamples_expr}"

def integrated_rowranges(wc=None):
    return ("resources/integrated_wgbs_hdf5/integrated_wgbs_combat_rowRanges.rds" if BC_METH
            else "resources/integrated_wgbs_hdf5/integrated_wgbs_rowRanges.rds")

rule h5_cpg_positions:
    output:
        "resources/integrated_wgbs_hdf5/h5_cpg_positions" + ("_combat" if BC_METH else "") + ".tsv"
    input:
        integrated_rowranges
    conda:
        "../envs/dataset_integration.yml"
    log:
        "logs/stability_selection/h5_cpg_positions.log"
    shell:
        "Rscript workflow/scripts/h5_cpg_positions.R {input} {output} > {log} 2>&1"

rule stability_kbest_features:   # cluster_cpgs: 1
    output:
        sel  = STAB_DIR + "/subsample_selections_shard{shard}.tsv",
        cand = STAB_DIR + "/candidates_shard{shard}.tsv"
    input:
        sel_check = "results/checks/model_selection" + TAG + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}",
        coldata   = "resources/integrated_wgbs_hdf5/meth_colData.tsv"
    params:
        sel_dir  = "results/sel_feats" + TAG,
        xgb_dir  = "data/xgb_clustered" + TAG,
        prefix   = stab_prefix,
        out_dir  = STAB_DIR,
        n_shards = STAB_SHARDS,
        n_pairs  = STAB_PAIRS
    conda:
        "../envs/py_combat.yml"
    log:
        "logs/stability_selection/kbest_features_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}_shard{shard}.log"
    shell:
        """
        python workflow/scripts/stability_selection_features.py \
        --sel_dir={params.sel_dir} --xgb_dir={params.xgb_dir} --prefix={params.prefix} \
        --coldata={input.coldata} --n_pairs={params.n_pairs} \
        --shard={wildcards.shard}/{params.n_shards} --out_dir={params.out_dir} > {log} 2>&1
        """

rule stability_kbest:   # cluster_cpgs: 0
    output:
        sel  = STAB_DIR + "/subsample_selections_shard{shard}.tsv",
        cand = STAB_DIR + "/candidates_shard{shard}.tsv"
    input:
        sel_check = "results/checks/model_selection" + TAG + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}",
        h5        = integrated_h5,
        positions = rules.h5_cpg_positions.output,
        coldata   = "resources/integrated_wgbs_hdf5/meth_colData.tsv",
        counts    = counts_train
    params:
        sel_dir  = "results/sel_feats" + TAG,
        cand_dir = "resources/CpG_data/top_K",
        prefix   = stab_prefix,
        out_dir  = STAB_DIR,
        n_shards = STAB_SHARDS,
        n_pairs  = STAB_PAIRS
    conda:
        "../envs/py_combat.yml"
    log:
        "logs/stability_selection/kbest_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}_shard{shard}.log"
    shell:
        """
        python analyses/scripts/stability_selection_kbest.py \
        --sel_dir={params.sel_dir} --cand_dir={params.cand_dir} --h5={input.h5} \
        --positions={input.positions} --coldata={input.coldata} --counts_train={input.counts} \
        --prefix={params.prefix} --n_pairs={params.n_pairs} \
        --shard={wildcards.shard}/{params.n_shards} --out_dir={params.out_dir} > {log} 2>&1
        """

rule stability_summary:
    output:
        genes = STAB_DIR + "/stability_genes.tsv",
        freq  = STAB_DIR + "/cpg_freq.tsv"
    input:
        sel  = expand(STAB_DIR + "/subsample_selections_shard{shard}.tsv", shard=range(STAB_SHARDS), allow_missing=True),
        cand = expand(STAB_DIR + "/candidates_shard{shard}.tsv", shard=range(STAB_SHARDS), allow_missing=True)
    params:
        in_dir = STAB_DIR
    conda:
        "../envs/r_stability.yml"
    log:
        "logs/stability_selection/summary_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.log"
    shell:
        "Rscript analyses/scripts/stability_selection_summary.R {params.in_dir} " + str(STAB_PAIRS) + " " + str(STAB_PFER) + " > {log} 2>&1"

rule stability_regions_corr:
    output:
        STAB_DIR + "/region_map_corr" + STAB_CORR + "_shard{shard}.tsv"
    input:
        cand      = STAB_DIR + "/candidates_shard{shard}.tsv",
        h5        = integrated_h5,
        positions = rules.h5_cpg_positions.output,
        coldata   = "resources/integrated_wgbs_hdf5/meth_colData.tsv",
        train     = "resources/split/train_samples.txt"
    params:
        in_dir   = STAB_DIR,
        cand_dir = "resources/CpG_data/top_K",
        prefix   = stab_prefix
    conda:
        "../envs/py_combat.yml"
    log:
        "logs/stability_selection/corr_regions_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}_shard{shard}.log"
    shell:
        """
        python analyses/scripts/stability_regions_corr.py \
        --in_dir={params.in_dir} --cand_dir={params.cand_dir} --h5={input.h5} \
        --positions={input.positions} --coldata={input.coldata} --train_samples={input.train} \
        --thresholds=""" + STAB_CORR + " --max_gap=" + str(STAB_GAP) + """ --prefix={params.prefix} \
        --shard={wildcards.shard} > {log} 2>&1
        """

rule stability_regions:
    output:
        freq  = STAB_DIR + "/region_freq_corr" + STAB_CORR + ".tsv",
        genes = STAB_DIR + "/region_genes_corr" + STAB_CORR + ".tsv"
    input:
        maps = expand(STAB_DIR + "/region_map_corr" + STAB_CORR + "_shard{shard}.tsv", shard=range(STAB_SHARDS), allow_missing=True),
        sel  = expand(STAB_DIR + "/subsample_selections_shard{shard}.tsv", shard=range(STAB_SHARDS), allow_missing=True),
        cand = expand(STAB_DIR + "/candidates_shard{shard}.tsv", shard=range(STAB_SHARDS), allow_missing=True)
    params:
        in_dir   = STAB_DIR,
        cand_dir = "resources/CpG_data/top_K"
    threads: 12
    conda:
        "../envs/r_stability.yml"
    log:
        "logs/stability_selection/regions_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.log"
    shell:
        "Rscript analyses/scripts/stability_selection_regions.R {params.in_dir} {params.cand_dir} "
        + str(STAB_GAP) + " " + str(STAB_PAIRS) + " " + str(STAB_PFER) + " {threads} corr" + STAB_CORR + " > {log} 2>&1"

if CLUSTER:
    ruleorder: stability_kbest_features > stability_kbest
else:
    ruleorder: stability_kbest > stability_kbest_features

STABLE_SEL_DIR = ("results/sel_feats" + TAG + ("_stable" if CLUSTER else "_stable_corr" + STAB_CORR)
                  + "/{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}")

rule stable_sel_feats:
    output:
        directory(STABLE_SEL_DIR)
    input:
        freq    = STAB_DIR + "/cpg_freq.tsv",
        regions = [] if CLUSTER else STAB_DIR + "/region_freq_corr" + STAB_CORR + ".tsv"
    params:
        stab_dir = STAB_DIR,
        prefix   = lambda wc: stab_prefix(wc) + f"_{wc.CV}",
        cmd      = ("workflow/scripts/stable_sel_feats_features.py" if CLUSTER
                    else "analyses/scripts/stable_sel_feats.py --tag=corr" + STAB_CORR)
    conda:
        "../envs/py_combat.yml"
    log:
        "logs/stability_selection/stable_sel_feats_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.log"
    shell:
        """
        python {params.cmd} --stab_dir={params.stab_dir} --out_dir={output} --prefix={params.prefix} > {log} 2>&1
        """
