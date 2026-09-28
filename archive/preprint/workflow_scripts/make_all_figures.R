# ============================================================
# make_all_figures.R
# ------------------------------------------------------------
# Regenerates the publication figures for the Culex pipiens
# complex pangenome paper from results/ on disk.
#
# Runs two ways:
#   - as a Snakemake rule (paths come from snakemake@params)
#   - standalone: Rscript make_all_figures.R [project_dir]
#     or set CXP_PROJ, or edit proj_default below.
#
# Repo: https://github.com/tylermaire/cx_pipiens_pangenome
# ============================================================

# --- 1. Packages ----------------------------------------------------------
cran_pkgs <- c("data.table", "ggplot2", "dplyr", "tidyr", "patchwork",
               "scales", "jsonlite", "ggrepel", "ggforce", "gridExtra",
               "RColorBrewer")
missing <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
bioc_pkgs <- c("ggtree", "treeio")
missing_bioc <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

suppressPackageStartupMessages({
  invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))
})

# --- 2. Paths -------------------------------------------------------------
# Edit this if running standalone from somewhere other than the project root.
proj_default <- "."

proj <- if (exists("snakemake")) {
  snakemake@params[["proj"]]
} else if (nzchar(Sys.getenv("CXP_PROJ"))) {
  Sys.getenv("CXP_PROJ")
} else if (length(commandArgs(trailingOnly = TRUE))) {
  commandArgs(trailingOnly = TRUE)[1]
} else {
  proj_default
}

res_dir   <- file.path(proj, "results")
out_dir   <- file.path(proj, "figures")
# anvi'o is run separately in the anvi'o web interface, not by this workflow;
# drop its exported gene-cluster summary here to regenerate that figure.
anvio_dir <- file.path(proj, "anvio_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("project: ", normalizePath(proj, mustWork = FALSE))

samples_ingroup <- c("Cx_molestus", "Cx_pallens", "Cx_pipiens", "Cx_quinquefasciatus")
samples_all     <- c(samples_ingroup, "Cx_tarsalis")
species_pretty <- c(Cx_molestus = "Cx. molestus", Cx_pallens = "Cx. pallens",
                    Cx_pipiens = "Cx. pipiens", Cx_quinquefasciatus = "Cx. quinquefasciatus",
                    Cx_tarsalis = "Cx. tarsalis")

# --- 3. Palette (MDPI-friendly, colour-blind safe) -------------------------
pal_grey <- list(dark = "#444444", mid = "#999999", light = "#DDDDDD")

species_colors <- c(Cx_molestus = "#E69F00", Cx_pallens = "#56B4E9",
                    Cx_pipiens = "#009E73", Cx_quinquefasciatus = "#CC79A7",
                    Cx_tarsalis = "#999999")
comp_colors  <- c(core = "#08519C", shell = "#6BAED6", cloud = "#C6DBEF")
busco_colors <- c("Single copy" = "#1A9850", "Multi copy" = "#A6D96A",
                  "Fragmented" = "#FDAE61", "Missing" = "#D73027")
te_colors    <- c("Interspersed TE" = "#66C2A5", "Simple repeats" = "#FC8D62",
                  "Low complexity" = "#8DA0CB", "Unmasked" = "#B3B3B3")
chrom_colors <- c(chr1 = "#377EB8", chr2 = "#4DAF4A", chr3 = "#E41A1C")
pal <- list(dark_grey = pal_grey$dark, light_grey = pal_grey$light)

theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(panel.grid.minor = element_blank(),
          axis.title  = element_text(size = base_size - 1),
          axis.text   = element_text(size = base_size - 2, color = "black"),
          plot.title  = element_text(size = base_size, face = "plain"),
          legend.title= element_text(size = base_size - 1, face = "bold"),
          legend.text = element_text(size = base_size - 1))
}

save_fig <- function(p, name, width, height, dpi = 300) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height)
  message("  wrote ", name, ".{png,pdf}")
}

# Screen rendering is cosmetic and some devices choke on large patchworks;
# the files are written by save_fig regardless.
show_plot <- function(p) invisible(try(print(p), silent = TRUE))

# =========================================================================
# Figure 1 — Assembly quality (BUSCO proteins + QUAST stats table)
# =========================================================================
fig1_assembly_quality <- function() {
  message("Figure 1: Assembly quality")
  busco <- lapply(samples_all, function(s) {
    fp <- file.path(res_dir, "busco_proteins", s,
                    paste0("short_summary.specific.diptera_odb10.", s, ".json"))
    if (!file.exists(fp)) return(NULL)
    r <- jsonlite::fromJSON(fp)$results
    data.frame(sample = s,
               "Single copy" = r[["Single copy percentage"]],
               "Multi copy"  = r[["Multi copy percentage"]],
               "Fragmented"  = r[["Fragmented percentage"]],
               "Missing"     = r[["Missing percentage"]],
               check.names = FALSE)
  })
  busco_long <- pivot_longer(do.call(rbind, busco), -sample,
                             names_to = "category", values_to = "pct")
  busco_long$category <- factor(busco_long$category,
                                levels = c("Single copy","Multi copy","Fragmented","Missing"))
  busco_long$sample <- factor(busco_long$sample, levels = rev(samples_all))

  p1a <- ggplot(busco_long, aes(x = pct, y = sample, fill = category)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
    geom_text(data = busco_long[busco_long$pct >= 4, ],
              aes(label = sprintf("%.1f", pct)),
              position = position_stack(vjust = 0.5), size = 3,
              color = ifelse(busco_long$category[busco_long$pct >= 4] == "Single copy",
                             "white", "black")) +
    scale_fill_manual(values = busco_colors) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 100.5),
                       breaks = c(0, 25, 50, 75, 100)) +
    scale_y_discrete(labels = function(x) parse(text = paste0("italic('", species_pretty[x], "')"))) +
    labs(title = "(a) Protein-mode BUSCO completeness",
         x = expression("BUSCO (% of 3,285 " * italic("diptera_odb10") * " markers)"),
         y = NULL, fill = NULL) +
    theme_pub() +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank())

  quast <- fread(file.path(res_dir, "quast", "transposed_report.tsv"))
  quast$Assembly <- factor(quast$Assembly, levels = samples_all)
  tab <- quast[order(Assembly), .(
    Species = species_pretty[as.character(Assembly)],
    "Size (Mb)" = round(`Total length (>= 0 bp)` / 1e6, 1),
    "Contigs"  = `# contigs`,
    "N50 (Mb)" = round(N50 / 1e6, 1),
    "GC (%)"   = round(`GC (%)`, 2))]
  tab_grob <- gridExtra::tableGrob(tab, rows = NULL,
    theme = gridExtra::ttheme_minimal(
      core    = list(fg_params = list(cex = 0.85,
                     fontface = ifelse(grepl("^Cx", tab$Species), 3, 1))),
      colhead = list(fg_params = list(cex = 0.9, fontface = "bold"))))
  p1b <- wrap_elements(full = tab_grob) +
    labs(title = "(b) Assembly statistics") + theme_pub() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())

  out <- p1a + p1b + plot_layout(widths = c(1.3, 1))
  show_plot(out)
  save_fig(out, "Figure_1_assembly_quality", width = 13, height = 5)
}

# =========================================================================
# Figure 2 — Pangenome composition
# =========================================================================
fig2_pangenome <- function() {
  message("Figure 2: Pangenome composition")
  summ  <- fread(file.path(res_dir, "pangenome", "pangenome_summary.tsv"))
  parts <- fread(file.path(res_dir, "pangenome", "partitioned_orthogroups.tsv"))
  parts[, comp_main := sub("_Cx_.*_specific$", "", compartment)]

  main3 <- summ[compartment %in% c("core","shell","cloud")]
  main3 <- main3[match(c("core","shell","cloud"), compartment)]
  main3[, pct := n_orthogroups / sum(n_orthogroups) * 100]
  main3[, y_mid := cumsum(n_orthogroups) - n_orthogroups / 2]
  main3[, label_inside := ifelse(pct >= 10,
    paste0(tools::toTitleCase(compartment), "\n",
           format(n_orthogroups, big.mark=","), "\n(", sprintf("%.1f", pct), "%)"), "")]
  main3[, label_outside := ifelse(pct < 10,
    paste0(tools::toTitleCase(compartment), " ",
           format(n_orthogroups, big.mark=","), " (", sprintf("%.1f", pct), "%)"), "")]
  main3$compartment <- factor(main3$compartment, levels = c("cloud","shell","core"))

  p2a <- ggplot(main3, aes(x = 2, y = n_orthogroups, fill = compartment)) +
    geom_bar(stat = "identity", color = "white", width = 1) +
    coord_polar(theta = "y", start = 0) + xlim(0.4, 3.0) +
    scale_fill_manual(values = comp_colors, breaks = c("core","shell","cloud"),
                      labels = function(x) tools::toTitleCase(x)) +
    geom_text(aes(x = 2.0, y = y_mid, label = label_inside), size = 3.1,
              lineheight = 0.95, color = "white", fontface = "bold", show.legend = FALSE) +
    geom_text(aes(x = 2.75, y = y_mid, label = label_outside), size = 3.0,
              color = "black", show.legend = FALSE) +
    annotate("text", x = 0.4, y = 0,
             label = format(sum(main3$n_orthogroups), big.mark=","),
             size = 6, fontface = "bold") +
    annotate("text", x = 0.4, y = 0, label = "orthogroups", size = 3.0, vjust = 3) +
    labs(title = "(a) Pangenome partition", fill = NULL) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(size = 11, hjust = 0.5, face = "plain"),
          legend.position = "none", plot.margin = margin(2, 2, 2, 2))

  per_form <- rbindlist(lapply(samples_ingroup, function(s)
    data.table(species = s,
               core  = sum(parts[comp_main == "core",  ..s][[1]]),
               shell = sum(parts[comp_main == "shell", ..s][[1]]),
               cloud = sum(parts[comp_main == "cloud", ..s][[1]]))))
  pf_long <- melt(per_form, id.vars = "species",
                  variable.name = "compartment", value.name = "genes")
  pf_long$compartment <- factor(pf_long$compartment, levels = c("cloud","shell","core"))
  totals <- pf_long[, .(total = sum(genes)), by = species]

  p2b <- ggplot(pf_long, aes(x = species, y = genes, fill = compartment)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.4) +
    geom_text(data = totals, aes(x = species, y = total,
              label = format(total, big.mark=",")),
              vjust = -0.4, inherit.aes = FALSE, size = 3.2) +
    scale_fill_manual(values = comp_colors, breaks = c("core","shell","cloud"),
                      labels = function(x) tools::toTitleCase(x)) +
    scale_x_discrete(labels = function(x) parse(text = paste0("italic('", species_pretty[x], "')"))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = "(b) Genes per species by compartment", x = NULL, y = "Genes", fill = NULL) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")

  cs <- summ[grepl("_specific$", compartment)]
  cs[, species := sub("cloud_(Cx_.*)_specific", "\\1", compartment)]
  cs$species <- factor(cs$species, levels = samples_ingroup)

  p2c <- ggplot(cs, aes(x = species, y = n_orthogroups, fill = species)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.4, show.legend = FALSE) +
    geom_text(aes(label = n_orthogroups), vjust = -0.4, size = 3.5) +
    scale_fill_manual(values = species_colors) +
    scale_x_discrete(labels = function(x) parse(text = paste0("italic('", species_pretty[x], "')"))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "(c) Form-specific cloud orthogroups", x = NULL, y = "Orthogroups") +
    theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))

  out <- p2a + p2b + p2c + plot_layout(ncol = 3, widths = c(1.3, 1.1, 1.0))
  show_plot(out)
  save_fig(out, "Figure_2_pangenome", width = 14.5, height = 5.2)
}

# =========================================================================
# Figure 3 — Species tree with concordance
# =========================================================================
fig3_species_tree <- function() {
  message("Figure 3: Species tree")
  tr <- ape::read.tree(file.path(res_dir, "phylo", "concord.cf.tree"))
  tr$tip.label <- species_pretty[tr$tip.label]

  p <- ggtree(tr, size = 0.9) +
    geom_tiplab(fontface = "italic", offset = 0.002, size = 4.2) +
    geom_nodelab(aes(label = label), hjust = -0.10, vjust = -0.6,
                 fontface = "bold", size = 3.6, color = "#222222") +
    geom_treescale(x = 0, y = -0.4, width = 0.05, fontsize = 3.2,
                   linesize = 0.7, offset = 0.1) +
    xlim(0, max(ggtree(tr)$data$x) * 1.55) +
    labs(caption = paste0(
      "Unrooted four-taxon topology: the data support only the split\n",
      "(Cx. pallens, Cx. quinquefasciatus) | (Cx. pipiens, Cx. molestus).\n",
      "The displayed root is arbitrary; Cx. tarsalis was not in the SCO set.\n",
      "Node label: ultrafast bootstrap / gene concordance / site concordance.")) +
    theme(plot.caption = element_text(hjust = 0.5, size = 9,
                                      face = "italic", color = "#444444"),
          plot.margin = margin(15, 30, 15, 15))

  show_plot(p)
  save_fig(p, "Figure_3_species_tree", width = 9.5, height = 6)
}

# =========================================================================
# Figure 4 — Pairwise whole-genome ANI  (was Figure 5)
# =========================================================================
fig4_ani <- function() {
  message("Figure 4: ANI heatmap")
  df <- fread(file.path(res_dir, "synteny", "ani_matrix.tsv"))
  setnames(df, "sample", "row_sample")
  df <- df[row_sample %in% samples_ingroup, c("row_sample", samples_ingroup), with = FALSE]
  long <- melt(df, id.vars = "row_sample", variable.name = "col_sample", value.name = "ani")
  long$ani <- suppressWarnings(as.numeric(long$ani))
  long$row_sample <- factor(long$row_sample, levels = samples_ingroup)
  long$col_sample <- factor(long$col_sample, levels = samples_ingroup)
  long$label <- ifelse(long$row_sample == long$col_sample, "100",
                       ifelse(is.na(long$ani), "NA", sprintf("%.2f", long$ani)))
  long$cell_fill <- ifelse(long$row_sample == long$col_sample, NA, long$ani)

  p <- ggplot(long, aes(x = col_sample, y = row_sample)) +
    geom_tile(aes(fill = cell_fill), color = "white", linewidth = 0.6) +
    geom_tile(data = long[long$row_sample == long$col_sample, ],
              fill = "#e0e0e0", color = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 4, fontface = "bold",
              color = ifelse(long$row_sample == long$col_sample, "#444444", "black")) +
    scale_fill_distiller(palette = "YlGnBu", direction = -1, limits = c(92.5, 95.0),
                         name = "skani\nANI (%)", na.value = "#f0f0f0") +
    scale_x_discrete(labels = function(x) parse(text = paste0("italic('", species_pretty[x], "')"))) +
    scale_y_discrete(labels = function(x) parse(text = paste0("italic('", species_pretty[x], "')")),
                     limits = rev) +
    coord_equal() + labs(x = NULL, y = NULL) + theme_pub() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())

  show_plot(p)
  save_fig(p, "Figure_4_ani", width = 7, height = 5.5)
}

# =========================================================================
# Figure 5 — Pairwise synteny dotplots  (was Figure 6)
# =========================================================================
fig5_synteny <- function() {
  message("Figure 5: Synteny")
  count_header_lines <- function(fp) {
    con <- file(fp, "r"); on.exit(close(con))
    n <- 0
    repeat {
      ln <- readLines(con, n = 1, warn = FALSE)
      if (length(ln) == 0 || !startsWith(ln, "#")) break
      n <- n + 1
    }
    n
  }
  read_genes <- function(sample) {
    fp <- file.path(res_dir, "annotation", paste0(sample, "_liftoff.gff3"))
    d <- fread(fp, sep = "\t", header = FALSE, fill = TRUE,
               skip = count_header_lines(fp),
               col.names = c("chrom","src","kind","start","end","score","strand","phase","attrs"))
    d <- d[kind == "gene"]
    d[, gene_id := sub("^ID=([^;]+).*", "\\1", attrs)]
    d[, gene_id := sub("^gene-", "", gene_id)]
    d[, gene_id := sub("_[0-9]+$", "", gene_id)]
    d[, mid := (start + end) / 2]
    d[, sample := sample]
    d[, .(sample, gene_id, chrom, mid)]
  }

  anchors <- rbindlist(lapply(samples_ingroup, read_genes))
  w_chrom <- dcast(anchors, gene_id ~ sample, value.var = "chrom",
                   fun.aggregate = function(x) x[1])
  w_pos   <- dcast(anchors, gene_id ~ sample, value.var = "mid",
                   fun.aggregate = function(x) x[1])

  plot_pair <- function(sx, sy, panel_letter) {
    d <- data.table(x_chrom = w_chrom[[sx]], y_chrom = w_chrom[[sy]],
                    x_pos   = w_pos[[sx]],   y_pos   = w_pos[[sy]])
    d <- d[!is.na(x_chrom) & !is.na(y_chrom)]
    rank_x <- d[, .N, by = x_chrom][order(-N)][1:3]
    rank_y <- d[, .N, by = y_chrom][order(-N)][1:3]
    d <- d[x_chrom %in% rank_x$x_chrom & y_chrom %in% rank_y$y_chrom]
    x_map <- setNames(c("chr1","chr2","chr3"), rank_x$x_chrom)
    d[, x_label := factor(x_map[x_chrom], levels = c("chr1","chr2","chr3"))]
    # Orientation is normalised per chromosome pair, matching synteny_summary:
    # a descending diagonal is an assembly-orientation difference, not biology.
    for (lv in levels(d$x_label)) {
      idx <- which(d$x_label == lv)
      if (length(idx) > 2 && cor(d$x_pos[idx], d$y_pos[idx],
                                 method = "spearman") < 0) {
        d$y_pos[idx] <- max(d$y_pos[idx]) - d$y_pos[idx]
      }
    }
    ggplot(d, aes(x = x_pos/1e6, y = y_pos/1e6, color = x_label)) +
      geom_point(size = 0.35, alpha = 0.55) +
      scale_color_manual(values = chrom_colors, name = "Reference chromosome",
                         guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(title = paste0("(", panel_letter, ") ", species_pretty[sx],
                          " vs ", species_pretty[sy]),
           x = paste0(species_pretty[sx], " (Mb)"),
           y = paste0(species_pretty[sy], " (Mb)")) +
      theme_pub(base_size = 10) +
      theme(axis.title = element_text(face = "italic", size = 9))
  }

  pairs_list <- list(c("Cx_pallens","Cx_pipiens"), c("Cx_pallens","Cx_quinquefasciatus"),
                     c("Cx_molestus","Cx_pallens"), c("Cx_molestus","Cx_pipiens"),
                     c("Cx_pipiens","Cx_quinquefasciatus"),
                     c("Cx_molestus","Cx_quinquefasciatus"))
  plots <- mapply(function(p, l) plot_pair(p[1], p[2], l),
                  pairs_list, letters[seq_along(pairs_list)], SIMPLIFY = FALSE)
  out <- wrap_plots(plots, ncol = 3, nrow = 2) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  show_plot(out)
  save_fig(out, "Figure_5_synteny", width = 13, height = 8.5)
}

# =========================================================================
# Figure 6 — TE composition  (was Figure 7)
# =========================================================================
fig6_te_composition <- function() {
  message("Figure 6: TE composition")
  rows <- lapply(samples_ingroup, function(s) {
    tl <- readLines(file.path(res_dir, "repeats", s, paste0(s, ".fasta.tbl")))
    num <- function(pat, lines) as.numeric(sub(pat, "\\1", grep(pat, lines, value = TRUE)[1]))
    total <- as.numeric(sub(".*total length:\\s+(\\d+)\\s+bp.*", "\\1",
                            grep("^total length:", tl, value = TRUE)[1]))
    inter <- as.numeric(sub(".*Total interspersed repeats:\\s+(\\d+)\\s+bp.*", "\\1",
                            grep("Total interspersed repeats:", tl, value = TRUE)[1]))
    simple <- as.numeric(sub(".*Simple repeats:\\s+\\d+\\s+(\\d+)\\s+bp.*", "\\1",
                             grep("^Simple repeats:", tl, value = TRUE)[1]))
    low <- as.numeric(sub(".*Low complexity:\\s+\\d+\\s+(\\d+)\\s+bp.*", "\\1",
                          grep("^Low complexity:", tl, value = TRUE)[1]))
    data.table(sample = s,
               "Interspersed TE" = inter/total*100, "Simple repeats" = simple/total*100,
               "Low complexity" = low/total*100,
               "Unmasked" = (total-inter-simple-low)/total*100, check.names = FALSE)
  })
  long <- melt(rbindlist(rows), id.vars = "sample",
               variable.name = "category", value.name = "pct")
  long$category <- factor(long$category,
    levels = c("Interspersed TE","Simple repeats","Low complexity","Unmasked"))
  long$sample <- factor(long$sample, levels = rev(samples_ingroup))

  p <- ggplot(long, aes(x = pct, y = sample, fill = category)) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE),
             color = "white", linewidth = 0.3) +
    geom_text(data = long[long$pct > 4, ], aes(label = sprintf("%.1f%%", pct)),
              position = position_stack(vjust = 0.5, reverse = TRUE), size = 3.2,
              color = ifelse(long$category[long$pct > 4] == "Interspersed TE",
                             "white", "black")) +
    scale_fill_manual(values = te_colors,
      breaks = c("Interspersed TE","Simple repeats","Low complexity","Unmasked")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 100.5)) +
    scale_y_discrete(labels = function(x) parse(text = paste0("italic('", species_pretty[x], "')"))) +
    labs(x = "Genome fraction (%)", y = NULL, fill = NULL) +
    theme_pub() +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank())

  show_plot(p)
  save_fig(p, "Figure_6_te_composition", width = 10, height = 4.8)
}

# =========================================================================
# Figure 7 — Anvi'o pangenome  (was Figure 8)
# Produced in the anvi'o web interface, not by this workflow. Export the
# gene-cluster summary into anvio_data/ to regenerate it.
# =========================================================================
fig7_anvio <- function() {
  message("Figure 7: Anvi'o pangenome")
  fp <- file.path(anvio_dir, "Culex_pipiens_complex_gene_clusters_summary.txt")
  if (!file.exists(fp)) {
    message("  SKIPPED: ", fp, " not found (anvi'o is run separately).")
    return(invisible(NULL))
  }
  d <- fread(fp, select = c("gene_cluster_id", "genome_name"))
  pa <- dcast(d[, .N, by = .(gene_cluster_id, genome_name)],
              gene_cluster_id ~ genome_name, value.var = "N", fill = 0)
  for (s in samples_ingroup) if (!s %in% names(pa)) pa[[s]] <- 0
  pa <- pa[, c("gene_cluster_id", samples_ingroup), with = FALSE]
  for (s in samples_ingroup) pa[[s]] <- as.integer(pa[[s]] > 0)
  pa[, pattern := apply(.SD, 1, paste, collapse = ""), .SDcols = samples_ingroup]
  pat <- pa[, .N, by = pattern][order(-N)]
  pat[, x := .I]

  p_bars <- ggplot(pat, aes(x = x, y = N)) +
    geom_bar(stat = "identity", fill = pal$dark_grey, color = "white", linewidth = 0.3) +
    geom_text(aes(label = format(N, big.mark = ",")), vjust = -0.4, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "Gene clusters", title = "(a) Cluster intersections") +
    theme_pub() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank())

  mat <- rbindlist(lapply(seq_len(nrow(pat)), function(i)
    data.table(x = i, species = samples_ingroup,
               present = as.integer(strsplit(pat$pattern[i], "")[[1]]))))
  mat$species <- factor(mat$species, levels = rev(samples_ingroup))
  mat$y <- as.integer(mat$species)
  conns <- mat[present == 1, .(ymin = min(y), ymax = max(y)), by = x][ymin < ymax]

  p_mat <- ggplot() +
    geom_segment(data = conns, aes(x = x, xend = x, y = ymin, yend = ymax),
                 color = pal$dark_grey, linewidth = 0.8) +
    geom_point(data = mat[present == 1, ], aes(x = x, y = y, color = species), size = 3.6) +
    geom_point(data = mat[present == 0, ], aes(x = x, y = y),
               color = pal$light_grey, size = 3.6) +
    scale_color_manual(values = species_colors, guide = "none") +
    scale_y_continuous(breaks = seq_along(samples_ingroup),
      labels = function(b) parse(text = paste0("italic('",
                                 species_pretty[rev(samples_ingroup)[b]], "')"))) +
    scale_x_continuous(breaks = NULL) +
    labs(x = "Intersection (sorted by size, left = largest)", y = NULL) +
    theme_void(base_size = 10) +
    theme(axis.text.y = element_text(hjust = 1), axis.title.x = element_text(size = 10),
          plot.margin = margin(2, 5, 5, 5))

  per_g <- data.table(species = factor(samples_ingroup, levels = samples_ingroup),
                      n = sapply(samples_ingroup, function(s) sum(pa[[s]])))
  p_pg <- ggplot(per_g, aes(x = n, y = species, fill = species)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.3, show.legend = FALSE) +
    geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.15, size = 3.3) +
    scale_fill_manual(values = species_colors) +
    scale_y_discrete(labels = function(x) parse(text = paste0("italic('", species_pretty[x], "')")),
                     limits = rev) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = "Clusters per genome", y = NULL, title = "(b) Per-genome cluster count") +
    theme_pub()

  sizes <- d[, .N, by = gene_cluster_id]
  p_sz <- ggplot(sizes, aes(x = N)) +
    geom_histogram(bins = 30, fill = pal$dark_grey, color = "white", linewidth = 0.3) +
    scale_x_log10() + annotation_logticks(sides = "b", size = 0.4) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("median %d | max %d", median(sizes$N), max(sizes$N)),
             hjust = 1.05, vjust = 1.5, size = 3.3, color = "#444") +
    labs(x = "Genes per cluster", y = "# clusters", title = "(c) Cluster-size distribution") +
    theme_pub()

  out <- (p_bars / p_mat + plot_layout(heights = c(3.5, 1.5))) |
         (p_pg / p_sz + plot_layout(heights = c(1, 1)))
  out <- out + plot_layout(widths = c(2.3, 1))
  show_plot(out)
  save_fig(out, "Figure_7_anvio_pangenome", width = 14, height = 8.5)
}

# =========================================================================
# Figure 8 — Key gene families heatmap  (was Figure 9)
# =========================================================================
fig8_key_families <- function() {
  message("Figure 8: Key gene families")
  d <- fread(file.path(res_dir, "functional", "key_families_wide.tsv"))
  mat <- as.matrix(d[, ..samples_ingroup])
  family_order <- d$family[order(-rowMeans(mat))]
  long <- melt(d, id.vars = "family", variable.name = "sample", value.name = "n")
  long$family <- factor(long$family, levels = family_order)
  long$sample <- factor(long$sample, levels = samples_ingroup)

  p <- ggplot(long, aes(x = sample, y = family, fill = n)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = n), size = 4) +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Copy\nnumber") +
    scale_x_discrete(labels = function(x) parse(text = paste0("italic('", species_pretty[x], "')"))) +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = "Gene family") + theme_pub() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), panel.grid = element_blank())

  show_plot(p)
  save_fig(p, "Figure_8_key_families", width = 8.5, height = 5.5)
}

# =========================================================================
# Main
# =========================================================================
# CAFE per-lineage expansion/contraction (previously Figure 4) is not plotted:
# the per-lineage signal is confounded with Liftoff annotation transfer
# (see rule cafe_transfer_bias), so it is reported as a supplementary table
# without lineage-specific interpretation.
main <- function() {
  cat("==> Generating figures into:", out_dir, "\n\n")
  figs <- list(fig1_assembly_quality, fig2_pangenome, fig3_species_tree,
               fig4_ani, fig5_synteny, fig6_te_composition,
               fig7_anvio, fig8_key_families)
  failed <- character(0)
  for (i in seq_along(figs)) {
    ok <- tryCatch({ figs[[i]](); TRUE },
                   error = function(e) { message("  ERROR: ", conditionMessage(e)); FALSE })
    if (!ok) failed <- c(failed, paste0("fig", i))
  }
  if (length(failed)) {
    cat("\n==> Finished with failures in:", paste(failed, collapse = ", "), "\n")
  } else {
    cat("\n==> All figures written.\n")
  }
}

if (!interactive()) main()
