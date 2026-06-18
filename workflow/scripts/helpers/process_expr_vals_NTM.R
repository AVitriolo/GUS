source("workflow/scripts/helpers/donorm.R")
source("workflow/scripts/helpers/filterGenes.R")

process_expr_vals_NTM <- function(input_expr_values_path_NT,
                                   input_expr_values_path_M,
                                   minCount_expr,
                                   minSamples_expr) {

    # ---- load and clean ----
    load_and_clean <- function(path) {
        df <- read.table(path, head = TRUE)
        df$transcript_id <- gsub("\\.[0-9]*", "", df$transcript_id)
        rownames(df) <- df$transcript_id
        df$transcript_id <- NULL
        df
    }

    NT <- load_and_clean(input_expr_values_path_NT)
    M  <- load_and_clean(input_expr_values_path_M)

    # ---- split NT into N and T by column name prefix ----
    N <- NT[, grepl("^N", colnames(NT)), drop = FALSE]
    T <- NT[, grepl("^T", colnames(NT)), drop = FALSE]

    cat(sprintf("Samples - N: %d, T: %d, M: %d\n", ncol(N), ncol(T), ncol(M)))

    # ---- filter by library size ----
    N.f <- N[, colSums(N) >= 5000000, drop = FALSE]
    T.f <- T[, colSums(T) >= 5000000, drop = FALSE]
    M.f <- M[, colSums(M) >= 5000000, drop = FALSE]

    # ---- filter by expression per condition ----
    N.ff <- filterGenes(N.f, minCount = minCount_expr, minSamples = minSamples_expr)
    T.ff <- filterGenes(T.f, minCount = minCount_expr, minSamples = minSamples_expr)
    M.ff <- filterGenes(M.f, minCount = minCount_expr, minSamples = minSamples_expr)

    # ---- normalize per condition ----
    N.tmm <- donorm(N.ff)
    T.tmm <- donorm(T.ff)
    M.tmm <- donorm(M.ff)

    # ---- intersect TxIDs across all 3 conditions ----
    common_txids <- intersect(intersect(rownames(N.tmm), rownames(T.tmm)), rownames(M.tmm))
    cat(sprintf("TxIDs passing filter - N: %d, T: %d, M: %d, common: %d\n",
                nrow(N.tmm), nrow(T.tmm), nrow(M.tmm), length(common_txids)))

    # ---- subset to common and merge ----
    merged <- cbind(
        N.tmm[common_txids, , drop = FALSE],
        T.tmm[common_txids, , drop = FALSE],
        M.tmm[common_txids, , drop = FALSE]
    )

    return(merged)
}
