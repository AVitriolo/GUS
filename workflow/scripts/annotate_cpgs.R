#!/usr/bin/env Rscript
# ============================================================
# annotate_cpgs.R
# Annotate two CpG sets into genomic categories (Promoter / Enhancer /
# Exon / Intron / Intergenic) using ENCODE cCREs + GENCODE, and compare
# their distributions (counts, density per Mb, % within stage).
#
#   set A = all ML-selected CpGs        (--bed_all)
#   set B = MethylDriver-significant    (--bed_md)
#
# Priority when a CpG overlaps several: Promoter > Enhancer_prox >
#   Enhancer_dist > Exon > Intron > Intergenic.
#
# Usage:
#   Rscript annotate_cpgs.R \
#     --bed_all=data/to_methyldriver/selected_features_percpg_hg38_v38_pc_50000_25_10_3.bed \
#     --bed_md=data/to_methyldriver/methyldriver_sig_cpgs.bed \
#     --pls=resources/regulatory_encode/GRCh38-cCREs.PLS.bed \
#     --els=resources/regulatory_encode/GRCh38-cCREs.ELS.bed \
#     --gtf=resources/reference_data/gencode/gencode.v38.annotation.gtf \
#     --out_dir=diagnostics/annotation
# ============================================================
suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(dplyr); library(ggplot2)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
bed_all <- args$bed_all
bed_md  <- args$bed_md
pls     <- args$pls
els     <- args$els
gtf_f   <- args$gtf
out_dir <- args$out_dir
label_all <- if (!is.null(args$label_all)) args$label_all else "All ML-selected"
label_md  <- if (!is.null(args$label_md))  args$label_md  else "MethylDriver sig"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

level_order <- c("Promoter","Enhancer_proximal","Enhancer_distal","Exon","Intron","Intergenic")
priority <- c(Promoter=1, Enhancer_proximal=2, Enhancer_distal=3, Exon=4, Intron=5)

# ---- regulatory / gene features (loaded once) ----
message("loading ENCODE cCREs + GENCODE ...")
promoters_gr <- GRanges(read.table(pls, header=FALSE, sep="\t")[,1:3] %>%
                        setNames(c("chr","start","end")))
promoters_gr$feature <- "Promoter"

enh_raw <- read.table(els, header=FALSE, sep="\t")[,c(1,2,3,6)] %>%
  setNames(c("chr","start","end","feature")) %>%
  mutate(feature = dplyr::recode(feature, "dELS"="Enhancer_distal", "pELS"="Enhancer_proximal"))
enhancers_gr <- GRanges(seqnames=enh_raw$chr, ranges=IRanges(enh_raw$start, enh_raw$end))
enhancers_gr$feature <- enh_raw$feature

gtf <- import(gtf_f, format="GTF", feature.type=c("exon","transcript"))
exons_gr <- gtf[gtf$type=="exon"]; exons_gr$feature <- "Exon"
transcripts_gr <- gtf[gtf$type=="transcript"]
introns_gr <- GenomicRanges::setdiff(transcripts_gr, exons_gr)
introns_gr$feature <- "Intron"

# ---- bp per category (for density) ----
category_bp <- function(gr) sum(as.numeric(width(GenomicRanges::reduce(gr))))
all_annotated_gr <- GenomicRanges::reduce(c(granges(promoters_gr), granges(enhancers_gr),
                                            granges(exons_gr), granges(introns_gr)))
bp_table <- data.frame(
  annotation = level_order,
  total_bp = c(category_bp(promoters_gr),
               category_bp(enhancers_gr[enhancers_gr$feature=="Enhancer_proximal"]),
               category_bp(enhancers_gr[enhancers_gr$feature=="Enhancer_distal"]),
               category_bp(exons_gr), category_bp(introns_gr),
               3.1e9 - sum(as.numeric(width(all_annotated_gr)))))

# ---- annotate one CpG BED ----
get_hits <- function(cpg_gr, feature_gr) {
  ov <- findOverlaps(cpg_gr, feature_gr)
  data.frame(cpg_id = cpg_gr$cpg_id[queryHits(ov)],
             annotation = as.character(feature_gr$feature[subjectHits(ov)]),
             stringsAsFactors = FALSE)
}

annotate_bed <- function(bed_path) {
  cpgs <- read.table(bed_path, header=FALSE, sep="\t")
  cpgs <- cpgs[,1:3]; colnames(cpgs) <- c("chr","start","end")
  cpgs$cpg_id <- seq_len(nrow(cpgs))
  gr <- GRanges(seqnames=cpgs$chr, ranges=IRanges(cpgs$start, cpgs$end), cpg_id=cpgs$cpg_id)

  best <- bind_rows(get_hits(gr, promoters_gr), get_hits(gr, enhancers_gr),
                    get_hits(gr, exons_gr), get_hits(gr, introns_gr)) %>%
    mutate(priority = priority[annotation]) %>%
    group_by(cpg_id) %>% slice_min(priority, n=1, with_ties=FALSE) %>% ungroup() %>%
    dplyr::select(cpg_id, annotation)

  res <- cpgs %>% left_join(best, by="cpg_id") %>%
    mutate(annotation = ifelse(is.na(annotation), "Intergenic", annotation))

  res %>% count(annotation, name="n_cpgs") %>%
    left_join(bp_table, by="annotation") %>%
    mutate(cpgs_per_mb = n_cpgs / (total_bp/1e6))
}

message("annotating set A (all) ..."); ca <- annotate_bed(bed_all)
message("annotating set B (md) ...");  cb <- annotate_bed(bed_md)

# ---- combine, collapse Enhancer_proximal into the plot set as in original ----
prep <- function(df) {
  df %>% filter(annotation != "Enhancer_proximal") %>%
    mutate(annotation = dplyr::recode(annotation, "Enhancer_distal"="Distal enhancer") %>%
             factor(levels=c("Promoter","Distal enhancer","Exon","Intron","Intergenic")))
}
comb <- bind_rows(prep(ca) %>% mutate(stage=label_all),
                  prep(cb) %>% mutate(stage=label_md)) %>%
  mutate(stage = factor(stage, levels=c(label_all, label_md)))

# percentage within each stage
tot <- comb %>% group_by(stage) %>% summarise(t=sum(n_cpgs), td=sum(cpgs_per_mb), .groups="drop")
comb <- comb %>% left_join(tot, by="stage") %>%
  mutate(pct = 100*n_cpgs/t, pct_dens = 100*cpgs_per_mb/td)

write.csv(comb, file.path(out_dir, "annotation_summary.csv"), row.names=FALSE)

pal <- c("#A8C8E8", "#00507A"); names(pal) <- c(label_all, label_md)

mk <- function(y, ylab, ttl, lab_fmt) {
  ggplot(comb, aes(x=annotation, y=.data[[y]], fill=stage)) +
    geom_col(position=position_dodge(width=0.7), width=0.65) +
    geom_text(aes(label=lab_fmt(.data[[y]])), position=position_dodge(width=0.7),
              vjust=-0.4, size=3.2) +
    scale_fill_manual(values=pal, name=NULL) +
    scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
    labs(title=ttl, x=NULL, y=ylab) +
    theme_classic(base_size=13) +
    theme(axis.text.x=element_text(angle=30, hjust=1), legend.position="top")
}

pdf(file.path(out_dir, "annotation_counts.pdf"), width=7, height=5)
print(mk("n_cpgs", "Number of CpGs", "CpGs per genomic category", scales::comma)); dev.off()

pdf(file.path(out_dir, "annotation_density.pdf"), width=7, height=5)
print(mk("cpgs_per_mb", "CpGs per Mb", "CpG density per genomic category",
         function(x) round(x,1))); dev.off()

pdf(file.path(out_dir, "annotation_pct.pdf"), width=7, height=5)
print(mk("pct", "% of CpGs", "CpG distribution (% within stage)",
         function(x) paste0(round(x,1),"%"))); dev.off()

pdf(file.path(out_dir, "annotation_pct_density.pdf"), width=7, height=5)
print(mk("pct_dens", "% of density", "CpG density distribution (% within stage)",
         function(x) paste0(round(x,1),"%"))); dev.off()

cat("done ->", out_dir, "\n")
print(comb %>% dplyr::select(stage, annotation, n_cpgs, cpgs_per_mb, pct))