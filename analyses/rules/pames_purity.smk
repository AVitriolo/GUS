rule pames_purity:
    output:
        purity     = "resources/purity/pames_purity_NT.tsv",
        info_sites = "resources/purity/pames_informative_sites_NT.tsv"
    input:
        cpgea = config["cpgea_wgbs_location"][0]
    params:
        min_cov         = config["pames_min_cov"][0],
        min_sample_frac = config["pames_min_sample_frac"][0],
        min_mean_diff   = config["pames_min_mean_diff"][0],
        max_na_data     = config["pames_max_na_data"][0],
        percentile_low  = config["pames_percentile_low"][0],
        percentile_high = config["pames_percentile_high"][0],
        max_sites       = config["pames_max_sites"][0],
        min_distance    = config["pames_min_distance"][0],
        cores           = config["n_jobs_xgboost"][0]
    log:
        "logs/pames_purity/pames_purity.log"
    shell:
        """
        mamba run -n r_pames Rscript analyses/scripts/pames_purity.R \
        --cpgea={input.cpgea} \
        --min_cov={params.min_cov} \
        --min_sample_frac={params.min_sample_frac} \
        --min_mean_diff={params.min_mean_diff} \
        --max_na_data={params.max_na_data} \
        --percentile_low={params.percentile_low} \
        --percentile_high={params.percentile_high} \
        --max_sites={params.max_sites} \
        --min_distance={params.min_distance} \
        --cores={params.cores} \
        --output_purity={output.purity} \
        --output_info_sites={output.info_sites} 2> {log}
        """
