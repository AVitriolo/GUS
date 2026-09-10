rule methyldriver_neighbors:            # runs ONCE
    output:
        neighbors = "results/methyldriver/neighbors_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.rds",
        matrix    = "results/methyldriver/featmatrix_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.rds"
    input:
        bed = "data/to_methyldriver/selected_features_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.bed"
    params:
        bigwig_dir    = config["encode_bigwig_dir"],
        exclude_files = "ENCFF582VMO.bigWig,ENCFF582VMO.bigWig.1",
        clusters_to_keep = lambda wc: (f"--clusters_to_keep={config['clusters_to_keep']}" if config.get("clusters_to_keep") else ""),
        n = 100, n_pcs = 20
    threads: 8
    conda: "../envs/r_methyldriver.yml"
    log: "logs/methyldriver/neighbors_{assembly_code}_v{gencode_version}_{tss_subset}_{distance}_{minCount_expr}_{minSamples_expr}_{CV}.log"
    shell:
        """
        export TMPDIR=/home/Prostate_GAS/resources/beegfs_scratch/tmp_gimme
        mkdir -p "$TMPDIR"
        Rscript workflow/scripts/build_neighbors.R \
        --input_bed={input.bed} --bigwig_dir={params.bigwig_dir} \
        --exclude_files={params.exclude_files} {params.clusters_to_keep} \
        --output_neighbors={output.neighbors} --output_matrix={output.matrix} \
        --num_cores={threads} --n={params.n} --n_pcs={params.n_pcs} > {log} 2>&1
        """
