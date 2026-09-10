#!/usr/bin/env Rscript
set.seed(123); options(scipen = 999)
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(ggrepel) })
source("workflow/scripts/helpers/MethylDriver.R")

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
input_neighbors  <- args$input_neighbors
input_matrix     <- args$input_matrix
input_meth_diff  <- args$input_meth_diff
conversion_table <- args$conversion_table
output_results   <- args$output_results
output_volcano   <- args$output_volcano
output_full_table<- args$output_full_table
q_cutoff         <- if (!is.null(args$q)) as.numeric(args$q) else 0.05

for (f in c(output_results, output_volcano, output_full_table))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)

closest_regions_list <- readRDS(input_neighbors)
fm                   <- readRDS(input_matrix)

mean_meth_diff <- read.table(input_meth_diff, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
home_region_methyl_change <- setNames(mean_meth_diff$delta, mean_meth_diff$clusterID)

target_order <- names(closest_regions_list)
home_region_methyl_change <- home_region_methyl_change[target_order]

transcript_ids <- sub("_.*$", "", target_order)
conv <- read.table(conversion_table, header = TRUE, stringsAsFactors = FALSE)
regions_genes <- setNames(conv$gene_name[match(transcript_ids, conv$transcript_id)], target_order)

stopifnot(all(names(home_region_methyl_change) == target_order))

res <- MethylDriver(home_region_methyl_change = home_region_methyl_change,
                    regions_genes = regions_genes,
                    closest_regions_list = closest_regions_list,
                    q_value_cutoff = q_cutoff)

saveRDS(res, output_results)
write.table(res$full_results, output_full_table, sep = "\t", quote = FALSE, row.names = FALSE)

df <- res$full_results
df$significance <- ifelse(df$q_value < q_cutoff & df$home_region_methyl_change > 0, "Hyper",
                   ifelse(df$q_value < q_cutoff & df$home_region_methyl_change < 0, "Hypo", "NS"))
top_genes <- df[order(df$p_value), ][seq_len(min(55, nrow(df))), ]
pdf(output_volcano, width = 12, height = 12)
print(ggplot(df, aes(home_region_methyl_change, -log10(p_value), color = significance)) +
  geom_point(alpha = .7) +
  ggrepel::geom_text_repel(data = top_genes, aes(label = gene), color = "black", size = 6, max.overlaps = 100) +
  scale_color_manual(values = c(Hypo="steelblue3", Hyper="#F8766B", NS="grey")) +
  labs(x = "Methylation Change", y = "-log10(p-value)", title = "MethylDriver") +
  theme_minimal())
dev.off()