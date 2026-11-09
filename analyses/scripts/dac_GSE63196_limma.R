# GSE63196 (Agilent GPL16604, quantile-normalised log2): 22Rv1, DU-145, LNCaP, 3 control vs 3 DAC (2 ug/mL, 4 days).
# limma DAC vs control within each cell line; one probe per gene symbol (highest mean signal over all 18 arrays).
# Out: gene, line, logFC, t, P.Value, adj.P.Val, ctrl_mean (baseline log2 signal in that line's controls).
suppressMessages({library(limma); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
mat_f <- args[1]; gpl_f <- args[2]; out_f <- args[3]

hdr <- readLines(gzfile(mat_f))
titles <- scan(text = sub("^!Sample_title\t", "", grep("^!Sample_title", hdr, value = TRUE)), what = "", sep = "\t", quiet = TRUE)
x <- fread(cmd = paste("zcat", mat_f), skip = "ID_REF", fill = TRUE)
x <- x[ID_REF != "!series_matrix_table_end"]
X <- as.matrix(x[, -1]); rownames(X) <- x$ID_REF

g <- fread(cmd = paste("zcat", gpl_f, "| sed -n '/^!platform_table_begin/,/^!platform_table_end/p' | sed '1d;$d'"), quote = "")
g <- g[ControlType == 0 & symbol != "" & ID %in% rownames(X)]
X <- X[g$ID, ]
keep <- tapply(seq_len(nrow(X)), g$symbol, function(i) i[which.max(rowMeans(X[i, , drop = FALSE]))])
X <- X[unlist(keep), ]; rownames(X) <- names(keep)
cat(nrow(X), "genes after one-probe-per-symbol\n")

line <- sub(" .*", "", sub("^(Control|DAC treated) ", "", titles))
trt  <- factor(ifelse(grepl("^DAC", titles), "DAC", "Ctrl"), levels = c("Ctrl", "DAC"))
res <- rbindlist(lapply(unique(line), function(l) {
  i <- line == l
  fit <- eBayes(lmFit(X[, i], model.matrix(~trt[i])))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), line = l, logFC = tt$logFC, t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val,
             ctrl_mean = rowMeans(X[rownames(tt), i & trt == "Ctrl"]))
}))
print(res[, .(up_fdr05 = sum(adj.P.Val < 0.05 & logFC > 0), down_fdr05 = sum(adj.P.Val < 0.05 & logFC < 0)), by = line])
fwrite(res, out_f, sep = "\t")
