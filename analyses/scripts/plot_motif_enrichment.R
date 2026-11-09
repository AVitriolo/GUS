#!/usr/bin/env Rscript
# =============================================================================
# plot_motif_enrichment.R
# Top-N enriched known motifs (gimme_known_enrichment.py report) as a bar plot:
# bar length = enrichment at 1% FPR (bars start at 1 = no enrichment over the
# background), fill = significance (-log10 BH q; single-hue scale, capped).
# Motifs ranked by enrichment among q < 0.05, labelled by motif name only.
# --groups (nonredundant_motifs.py): only group representatives, motifs longer
# than the windows excluded.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

args   <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
report <- fread(args$report)
n_top  <- as.integer(if (is.null(args$n)) 50 else args$n)
out    <- args$output                       # .pdf; a .png and a .tsv are written next to it
q_cap  <- as.numeric(if (is.null(args$q_cap)) 200 else args$q_cap)   # colour scale cap (-log10 q)

pick <- function(x) unique(toupper(trimws(strsplit(x, ",", fixed = TRUE)[[1]])))
# genes of a set of motifs, most frequent first (the representative's own genes win ties)
rank_genes <- function(lists, own) {
    v <- unlist(lapply(lists, function(x) pick(x)))
    v <- v[nzchar(v)]
    if (!length(v)) return(character(0))
    tab <- table(v)
    names(tab)[order(-as.integer(tab), !(names(tab) %in% own), names(tab))]
}
d0 <- report[`q-value (BH)` < 0.05 & `Enr. at 1% FPR` > 1]
if (!is.null(args$groups)) {
    # one bar per group of similar motifs (nonredundant_motifs.py), labelled with the genes of
    # the whole group: direct annotations; indirect/predicted ones in brackets when the group has
    # no direct one; groups with no annotated gene are dropped
    G <- fread(args$groups)
    reps <- G[is_representative == TRUE & excluded_too_long == FALSE, .(Motif, motif_id, group_size)]
    d0 <- d0[reps, on = "Motif", nomatch = NULL]
    lab <- rbindlist(lapply(d0$motif_id, function(r) {
        m <- G[representative == r]
        own_d <- pick(fcoalesce(m[motif_id == r, factors_direct], ""))
        gd <- rank_genes(fcoalesce(m$factors_direct, ""), own_d)
        gi <- if (length(gd)) character(0) else rank_genes(fcoalesce(m$factors_indirect, ""), character(0))
        data.table(motif_id = r, genes = paste(head(if (length(gd)) gd else gi, 3), collapse = ", "), direct = length(gd) > 0)
    }))
    d0 <- d0[lab, on = "motif_id"]
} else {
    d0[, group_size := 1L]
    d0[, genes := mapply(function(a, b) paste(head(if (nzchar(a)) pick(a) else pick(b), 3), collapse = ", "),
                         fcoalesce(factors_direct, ""), fcoalesce(factors_indirect, ""))]
    d0[, direct := nzchar(fcoalesce(factors_direct, ""))]
}
d0 <- d0[nzchar(genes)]                                   # no annotated TF: dropped
d <- d0[order(-`Enr. at 1% FPR`, `q-value (BH)`)][seq_len(min(n_top, .N))]
d[, motif_name := sub("^GM\\.5\\.0\\.", "", sub("_[^_]*$", "", Motif))]
d[, name := fifelse(direct, genes, paste0("(", genes, ")"))]
d[, name := fifelse(duplicated(name) | duplicated(name, fromLast = TRUE), sprintf("%s [%s]", name, motif_name), name)]
d[, name := fifelse(group_size > 1, sprintf("%s  [+%d similar]", name, group_size - 1), name)]
d[, name := factor(name, levels = rev(name))]
d[, mlq := -log10(pmax(`q-value (BH)`, 1e-300))]

ink <- c(primary = "#0b0b0b", secondary = "#52514e", grid = "#e6e5e0")
pal <- c("#f3e1a9", "#e3c9a5", "#d1b0a0", "#b88f9a", "#a06e93", "#723287")   # light yellow -> purple
lim <- c(floor(min(d$mlq) / 10) * 10, q_cap)
brk <- pretty(lim, 4); brk <- brk[brk >= lim[1] & brk <= lim[2]]
lab <- ifelse(brk == q_cap, paste0("≥", q_cap), brk)

p <- ggplot(d, aes(y = name)) +
    geom_rect(aes(xmin = 1, xmax = `Enr. at 1% FPR`, ymin = as.integer(name) - 0.4, ymax = as.integer(name) + 0.4,
                  fill = pmin(mlq, q_cap))) +
    scale_fill_gradientn(colours = pal, limits = lim, breaks = brk, labels = lab, oob = scales::squish,
                         name = "-log10 q (BH)") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
    coord_cartesian(xlim = c(1, max(d$`Enr. at 1% FPR`) * 1.01)) +
    labs(x = if (is.null(args$xlab)) "Enrichment at 1% FPR (selected vs matched unselected CpG windows)" else args$xlab, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = ink[["grid"]], linewidth = 0.3),
          axis.text.y = element_text(size = 7.5, colour = ink[["primary"]]),
          axis.text.x = element_text(colour = ink[["secondary"]]),
          axis.title.x = element_text(colour = ink[["secondary"]], size = 9),
          legend.title = element_text(size = 9, colour = ink[["primary"]]),
          legend.text = element_text(size = 8, colour = ink[["secondary"]]),
          plot.background = element_rect(fill = "white", colour = NA))

h <- 1.0 + 0.17 * nrow(d)
ggsave(out, p, width = 7.5, height = h, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", out), p, width = 7.5, height = h, dpi = 250, bg = "white")
keep_cols <- intersect(c("Motif", "name", "genes", "direct", "group_size", "factors_direct", "factors_indirect", "% matches input", "% matches background",
                         "Enr. at 1% FPR", "ROC AUC", "q-value (BH)"), names(d))
fwrite(cbind(rank = seq_len(nrow(d)), d[, ..keep_cols]), sub("\\.pdf$", ".tsv", out), sep = "\t")
cat(sprintf("%d motifs -> %s\n", nrow(d), out))
