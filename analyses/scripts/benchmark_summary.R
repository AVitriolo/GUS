# Feature-selector benchmark, step 4: stability of each selector over the half-samples.
# For every gene x method x half-sample set (SKB/ENET: all 100 and the first 20; SUSIE: the first 20):
#   phi, phi_region      Nogueira stability index (stabm) at CpG and co-methylated-region level
#   n_stable_cpg/region  CpGs / regions with frequency >= the stabs cutoff for PFER <= 1
#                        (complementary pairs, unimodality, q = mean selection size rounded up)
#   pct_full_in_stable_region  % of the full-train selected CpGs lying in a stable region
# Writes stability_genes.tsv in <bench_dir>.
#
# Usage: Rscript benchmark_summary.R <bench_dir> [pfer=1] [cores=12]

suppressPackageStartupMessages({ library(stabm); library(stabs); library(data.table); library(parallel) })

args  <- commandArgs(trailingOnly = TRUE)
bdir  <- args[1]
pfer  <- if (length(args) >= 2) as.numeric(args[2]) else 1
cores <- if (length(args) >= 3) as.integer(args[3]) else 12L

rd <- function(pattern, names) rbindlist(lapply(list.files(bdir, pattern, full.names = TRUE), fread, sep = "\t",
                                                header = FALSE, col.names = names, colClasses = "character"))
inp <- rd("^bench_inputs_shard.*\\.tsv$", c("TxID", "what", "v"))[what %in% c("region")]
sel <- rd("^bench_selections.*_shard.*\\.tsv$", c("TxID", "method", "b", "idx"))
region <- lapply(split(inp$v, inp$TxID), function(s) as.integer(strsplit(s, ",")[[1]]))
sel <- split(sel, by = "TxID")

idx_of <- function(s) if (is.na(s) || s == "") integer(0) else as.integer(strsplit(s, ",")[[1]]) + 1L

stab <- function(sets, p, n_pairs) {
    q <- ceiling(mean(lengths(sets)))
    cutoff <- if (q < 1) NA_real_ else tryCatch(
        stabsel_parameters(p = p, q = q, PFER = pfer, B = n_pairs, sampling.type = "SS",
                           assumption = "unimodal")$cutoff, error = function(e) NA_real_)
    freq <- table(factor(unlist(sets), levels = seq_len(p))) / length(sets)
    phi  <- if (all(lengths(sets) == 0)) NA_real_ else stabilityNogueira(sets, p = p)
    list(phi = phi, cutoff = cutoff, q = q,
         stable = if (is.na(cutoff)) integer(0) else which(freq >= cutoff))
}

one_gene <- function(gid) {
    s   <- sel[[gid]]
    reg <- region[[gid]]
    p_cpg <- length(reg); p_reg <- max(reg)
    res <- list()
    for (m in unique(s$method)) {
        sm   <- s[method == m]
        full <- idx_of(sm[b == "full", idx])
        hv   <- sm[b != "full"][order(as.integer(b))]
        for (nb in unique(c(20L, nrow(hv)))) {
            if (nrow(hv) < nb || nb < 2) next
            sets <- lapply(hv$idx[seq_len(nb)], idx_of)
            cpg  <- stab(sets, p_cpg, nb / 2)
            rg   <- stab(lapply(sets, function(x) unique(reg[x])), p_reg, nb / 2)
            res[[length(res) + 1]] <- data.table(
                TxID = gid, method = m, n_halves = nb, phi = cpg$phi, phi_region = rg$phi,
                q_cpg = cpg$q, q_region = rg$q, cutoff_cpg = cpg$cutoff, cutoff_region = rg$cutoff,
                n_stable_cpg = length(cpg$stable), n_stable_region = length(rg$stable),
                n_full = length(full),
                pct_full_stable_cpg = if (length(full)) 100 * mean(full %in% cpg$stable) else NA_real_,
                pct_full_in_stable_region = if (length(full)) 100 * mean(reg[full] %in% rg$stable) else NA_real_)
        }
    }
    rbindlist(res)
}

out <- rbindlist(mclapply(names(sel), one_gene, mc.cores = cores))
fwrite(out, file.path(bdir, "stability_genes.tsv"), sep = "\t")
print(out[, .(genes = .N, phi = median(phi, na.rm = TRUE), phi_region = median(phi_region, na.rm = TRUE),
              n_stable_region = median(n_stable_region), pct_full_in_stable_region = median(pct_full_in_stable_region, na.rm = TRUE)),
          by = .(method, n_halves)])
