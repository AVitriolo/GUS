rule batch_correct_meth:
    output:
        h5        = "resources/integrated_wgbs_hdf5/integrated_wgbs_combat.h5",
        coldata   = "resources/integrated_wgbs_hdf5/integrated_wgbs_combat_colData.rds",
        rowranges = "resources/integrated_wgbs_hdf5/integrated_wgbs_combat_rowRanges.rds",
        params    = "resources/integrated_wgbs_hdf5/combat_params_meth.npz",
        plot_ds   = "results/plots/batch_correct/meth_combat_by_dataset.pdf",
        plot_cond = "results/plots/batch_correct/meth_combat_by_condition.pdf"
    input:
        h5        = "resources/integrated_wgbs_hdf5/integrated_wgbs.h5",
        coldata   = "resources/integrated_wgbs_hdf5/integrated_wgbs_colData.rds",
        rowranges = "resources/integrated_wgbs_hdf5/integrated_wgbs_rowRanges.rds",
        meta_cd   = "resources/integrated_wgbs_hdf5/meth_colData.tsv",
        meta_rr   = "resources/integrated_wgbs_hdf5/meth_rowRanges.tsv"
    conda:
        "../envs/py_combat.yml"
    log:
        "logs/batch_correct/meth_combat.log"
    shell:
        """
        PYTHONPATH=workflow/scripts/helpers python workflow/scripts/combat_correction_meth.py \
        --in_h5={input.h5} \
        --in_rowranges={input.meta_rr} \
        --in_coldata={input.meta_cd} \
        --out_h5={output.h5} \
        --out_params={output.params} \
        --out_plot_ds={output.plot_ds} \
        --out_plot_cond={output.plot_cond} 2> {log}
        cp {input.coldata}   {output.coldata}
        cp {input.rowranges} {output.rowranges}
        """
