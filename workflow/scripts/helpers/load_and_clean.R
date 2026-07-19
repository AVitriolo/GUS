load_and_clean <- function(path) {
    df <- read.table(path, header = TRUE)
    df$transcript_id <- gsub("\\.[0-9]*", "", df$transcript_id)
    rownames(df) <- df$transcript_id
    df$transcript_id <- NULL
    df
}
