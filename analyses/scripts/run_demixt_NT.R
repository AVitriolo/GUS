# =============================================================================
# run_demixt_NT.R
# Purified CPGEA Tumour expression with the DeMixT package (>= 2.0.0 from
# github.com/wwylab/DeMixT: 1.16.0's DeMixT_preprocessing() calls an undefined
# quantile_normalization(); Wang et al. iScience 2018), two-component model:
# Tumour = mix of tumour + normal prostate, CPGEA Normals as the reference.
#
#   per split (train, test), run separately so test never informs train:
#     reference N1 : TRAIN Normals (the same reference for both splits)
#     mixed     Y  : the split's Tumours
#     1. DeMixT_preprocessing on raw counts -> ~9000 transcripts by expression
#        SD, upper-quartile scaled; used only to estimate the tumour
#        proportion pi (DeMixT_GS, --ngene genes, --nspikein spike-ins)
#     2. DeMixT_S2 with that pi -> purified tumour expression (decovExprT) for
#        every pipeline transcript (--stage all)
#
#   Counts get +1 before DeMixT_S2. S2 fits log2(counts), so it drops every
#   transcript with a zero in any sample; its own alternative
#   (filter.option = 2, +1 pseudocount) is hard-coded off, and adding 1 here
#   reproduces it, keeping all transcripts in both splits. DeMixT_GS already
#   replaces zeros with 1 internally.
#
#   --stage gs runs step 1 only (parameter grid, see run_demixt_grid_NT.sh).
#   Outputs are tagged ns<nspikein>_ng<ngene>. Normals are not deconvolved.
# =============================================================================
suppressPackageStartupMessages({ library(DeMixT); library(SummarizedExperiment); library(data.table) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

counts_path  <- args$counts               # resources/RNA/kallisto_counts_NT.tsv (raw)
tx_path      <- args$transcripts          # pipeline transcript universe (first column)
coldata_path <- args$coldata              # integrated colData: condition / split
outdir       <- args$outdir
nthread      <- as.integer(if (is.null(args$nthread)) 20 else args$nthread)
ngene        <- as.integer(if (is.null(args$ngene)) 1500 else args$ngene)
nspikein     <- if (is.null(args$nspikein)) NULL else as.integer(args$nspikein)
stage        <- if (is.null(args$stage)) "all" else args$stage
splits       <- if (is.null(args$splits)) c("train", "test") else strsplit(args$splits, ",")[[1]]
given_pi     <- args$given_pi                 # optional: purity table (e.g. PAMES); skips DeMixT_GS
stopifnot(stage %in% c("gs", "all"))

set.seed(1234)   # DeMixT uses random initial values and spike-ins
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cd <- as.data.frame(readRDS(coldata_path))
tx <- fread(tx_path, select = 1)[[1]]
raw <- fread(counts_path)
raw <- raw[transcript_id %in% tx]
cm  <- as.matrix(raw[, -1]); rownames(cm) <- raw$transcript_id

normals_train <- rownames(cd)[cd$condition == "Normal" & cd$split == "train"]

for (s in splits) {
    t0      <- Sys.time()
    # "testjoint": fit ALL Tumours (train + test) together, keep only the test
    # Tumours' outputs. Test borrows the train Tumours for a stable fit; train
    # outputs still come from the train-only run, so no test information
    # reaches the training data.
    tumours  <- if (s == "testjoint") rownames(cd)[cd$condition == "Tumour"] else
                    rownames(cd)[cd$condition == "Tumour" & cd$split == s]
    keep_out <- if (s == "testjoint") rownames(cd)[cd$condition == "Tumour" & cd$split == "test"] else tumours
    ns       <- if (is.null(nspikein)) min(200, ceiling(length(tumours) * 0.3)) else nspikein
    tag     <- if (is.null(given_pi)) sprintf("%s_ns%d_ng%d", s, ns, ngene) else sprintf("%s_givenpi", s)
    raw_s   <- cm[, c(normals_train, tumours)]
    raw_s   <- raw_s[rowSums(raw_s) > 0, ]

    if (!is.null(given_pi)) {
        # purity supplied (same-block DNA/RNA in CPGEA): no proportion
        # estimation, DeMixT_S2 purifies with the given tumour share
        gp <- read.delim(given_pi, row.names = 1)
        gp <- setNames(gp[[grep("^purity$", colnames(gp), ignore.case = TRUE)]], rownames(gp))
        if (length(setdiff(tumours, names(gp)))) stop("No given purity for: ", paste(setdiff(tumours, names(gp)), collapse = ", "))
        gs <- list(pi = rbind(PiN1 = 1 - gp[tumours], PiT = gp[tumours]))
        colnames(gs$pi) <- tumours
        cat(sprintf("[%s] %d Tumours vs %d train Normals | given purity from %s\n",
                    tag, length(tumours), length(normals_train), basename(given_pi)))
    } else {
    # DeMixT >= 2.0: preprocessing takes RAW counts, selects transcripts by SD,
    # then upper-quartile scales the selection itself
    pre <- DeMixT_preprocessing(raw_s, normals_train, tumours,
                                cutoff_normal_range = c(0.1, 1.0),
                                cutoff_tumor_range  = c(0, 2.5), cutoff_step = 0.1)
    sel <- pre$count.matrix
    cat(sprintf("[%s] %d Tumours vs %d train Normals | %d transcripts, %d kept for pi (sd cutoffs N %s, T %s)\n",
                tag, length(tumours), length(normals_train), nrow(raw_s), nrow(sel),
                paste(pre$sd_cutoff_normal, collapse = "-"), paste(pre$sd_cutoff_tumor, collapse = "-")))

    # step 1 is the slow part: cache it so a failure in step 2 does not redo it
    gs_rds <- file.path(outdir, sprintf("demixt_GS_%s.rds", tag))
    if (file.exists(gs_rds)) {
        gs <- readRDS(gs_rds)
        cat(sprintf("[%s] reusing cached proportions %s\n", tag, gs_rds))
    } else {
        gs <- DeMixT_GS(data.Y  = SummarizedExperiment(assays = list(counts = sel[, tumours])),
                        data.N1 = SummarizedExperiment(assays = list(counts = sel[, normals_train])),
                        niter = 50, nbin = 50, if.filter = TRUE,
                        filter.sd = pre$sd_cutoff_normal[2],
                        ngene.Profile.selected = ngene, ngene.selected.for.pi = ngene,
                        mean.diff.in.CM = 0.25, nspikein = ns, tol = 10^(-5), nthread = nthread)
        saveRDS(gs, gs_rds)
    }
    }
    stopifnot(identical(colnames(gs$pi), tumours))    # spike-ins are dropped by DeMixT
    pi_t <- gs$pi["PiT", keep_out]
    fwrite(data.table(sample = keep_out, pi_tumour = pi_t),
           file.path(outdir, sprintf("demixt_pi_%s.tsv", tag)), sep = "\t")
    cat(sprintf("[%s] pi estimated (%.1f min): tumour proportion median %.2f, range %.2f-%.2f\n",
                tag, as.numeric(difftime(Sys.time(), t0, units = "mins")),
                median(pi_t), min(pi_t), max(pi_t)))
    if (stage == "gs") next

    m  <- scale_normalization_75th_percentile(raw_s)  # all transcripts
    m1 <- m + 1                                        # = DeMixT_S2 filter.option 2
    s2 <- DeMixT_S2(data.Y  = SummarizedExperiment(assays = list(counts = m1[, tumours])),
                    data.N1 = SummarizedExperiment(assays = list(counts = m1[, normals_train])),
                    data.N2 = NULL,
                    givenpi = c(t(gs$pi["PiN1", , drop = FALSE])),
                    nbin = 50, nthread = nthread)
    expr_t <- s2$decovExprT
    stopifnot(!is.null(rownames(expr_t)), identical(colnames(expr_t), tumours))
    cat(sprintf("[%s] DeMixT_S2 deconvolved %d of %d transcripts (dropped: constant across Normals)\n",
                tag, nrow(expr_t), nrow(m1)))

    expr_t <- expr_t[, keep_out, drop = FALSE]
    fwrite(data.table(transcript_id = rownames(expr_t), expr_t),
           file.path(outdir, sprintf("demixt_purified_tumour_%s.tsv.gz", tag)), sep = "\t")
    fwrite(data.table(transcript_id = rownames(m1), m1[, keep_out]),
           file.path(outdir, sprintf("bulk_scaled_plus1_tumour_%s.tsv.gz", tag)), sep = "\t")
    cat(sprintf("[%s] done in %.1f min -> %s\n", tag,
                as.numeric(difftime(Sys.time(), t0, units = "mins")), outdir))
}
