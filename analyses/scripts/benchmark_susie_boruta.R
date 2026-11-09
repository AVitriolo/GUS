# Feature-selector benchmark, step 2: SuSiE and Boruta on the same genes, data and half-samples as
# benchmark_selectors.py (reads its bench_inputs / bench_halves shards).
#   SUSIE   susieR::susie(X, y, L = 10); selected = CpGs in the 95% credible sets (default purity filter)
#   BORUTA  Boruta::Boruta(X, y) (ranger importance, pValue 0.01, maxRuns 100); selected = Confirmed
# Both on the full train set ("full"); SuSiE also on the first susie_halves half-samples, Boruta on the first
# boruta_halves half-samples only (cost). No covariate adjustment.
#
# Usage: Rscript benchmark_susie_boruta.R <bench_dir> <shard> <h5> [boruta_halves=0] [threads=1] [susie_halves=20]

suppressPackageStartupMessages({ library(susieR); library(Boruta); library(rhdf5); library(data.table) })

args   <- commandArgs(trailingOnly = TRUE)
bdir   <- args[1]
shard  <- args[2]
h5     <- args[3]
n_bor_halves <- if (length(args) >= 4) as.integer(args[4]) else 0L
threads      <- if (length(args) >= 5) as.integer(args[5]) else 1L
n_sus_halves <- if (length(args) >= 6) as.integer(args[6]) else 20L
suppressPackageStartupMessages(library(Rfast))                         # faster credible-set purity

inp    <- fread(file.path(bdir, sprintf("bench_inputs_shard%s.tsv", shard)), header = FALSE, sep = "\t",
                col.names = c("TxID", "what", "v"), colClasses = "character")
halves <- fread(file.path(bdir, sprintf("bench_halves_shard%s.tsv", shard)), header = FALSE, sep = "\t",
                col.names = c("TxID", "b", "idx"), colClasses = "character")
if (!nrow(halves)) halves <- data.table(TxID = character(0), b = character(0), idx = character(0))   # null run
rows   <- scan(file.path(bdir, "bench_train_rows.txt"), quiet = TRUE) + 1L
num    <- function(s) as.numeric(strsplit(s, ",")[[1]])

susie_select <- function(X, y) {
    fit <- tryCatch(susie(X, y, L = 10), error = function(e) NULL)
    if (is.null(fit) || is.null(fit$sets$cs)) return(integer(0))
    sort(unique(unlist(fit$sets$cs)))
}
boruta_select <- function(X, y) {
    b <- Boruta(X, y, pValue = 0.01, maxRuns = 100, num.threads = threads)
    which(b$finalDecision == "Confirmed")
}
out <- file(file.path(bdir, sprintf("bench_selections_susie_boruta_shard%s.tsv", shard)), "w")
emit <- function(TxID, method, b, idx) {
    writeLines(sprintf("%s\t%s\t%s\t%s", TxID, method, b, paste(idx - 1L, collapse = ",")), out)
    flush(out)
}

for (gid in unique(inp$TxID)) {
    g    <- inp[TxID == gid]
    cols <- as.integer(num(g[what == "cols", v])) + 1L
    y    <- num(g[what == "y_train", v])
    X    <- h5read(h5, "M", index = list(min(cols):max(cols), NULL))      # CpGs x samples
    X    <- t(X[cols - min(cols) + 1L, rows, drop = FALSE])
    colnames(X) <- paste0("c", seq_len(ncol(X)))
    H    <- halves[TxID == gid]
    set.seed(as.integer(sub("ENST0*", "", gid)) %% .Machine$integer.max)

    t0 <- Sys.time()
    emit(gid, "SUSIE", "full", susie_select(X, y))
    for (i in seq_len(min(n_sus_halves, nrow(H)))) {
        h <- as.integer(num(H$idx[i])) + 1L
        emit(gid, "SUSIE", H$b[i], susie_select(X[h, ], y[h]))
    }
    t1 <- Sys.time()
    emit(gid, "BORUTA", "full", boruta_select(X, y))
    for (i in seq_len(min(n_bor_halves, nrow(H)))) {
        h <- as.integer(num(H$idx[i])) + 1L
        emit(gid, "BORUTA", H$b[i], boruta_select(X[h, ], y[h]))
    }
    cat(sprintf("shard %s %s: susie %.0fs, boruta %.0fs\n", shard, gid,
                as.numeric(difftime(t1, t0, units = "secs")), as.numeric(difftime(Sys.time(), t1, units = "secs"))))
}
close(out)
