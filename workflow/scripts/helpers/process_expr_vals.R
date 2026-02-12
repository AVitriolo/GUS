source("workflow/scripts/helpers/donorm.R")
source("workflow/scripts/helpers/filterGenes.R")

process_expr_vals <- function(input_expr_values_path, minCount_expr, minSamples_expr){
    ExprValues <- read.table(input_expr_values_path,head=T) # "RNA/kallisto_counts.tsv"
    rownames(ExprValues) <- ExprValues$transcript_id
    ExprValues$transcript_id <- NULL
    ExprValues.f <- ExprValues[, which(colSums(ExprValues) >= 5000000)]
    ExprValues.ff <- filterGenes(ExprValues.f,minCount=minCount_expr,minSamples=minSamples_expr) # minCounts = 20, minSamples = 10
    ExprValues.tmm <- donorm(ExprValues.ff)
    return(ExprValues.tmm)
}