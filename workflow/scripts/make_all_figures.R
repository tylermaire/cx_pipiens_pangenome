# ============================================================
# make_all_figures.R
# ------------------------------------------------------------
# Master script: regenerates ALL 9 publication figures for the
# Culex pipiens complex pangenome paper from results/ on disk.
#
# Run from RStudio or `Rscript make_all_figures.R`.
# Outputs PNG + PDF for each figure to figures/.
#
# Author: Tyler Maire / paper companion repo
# Repo:   https://github.com/tylermaire/cx_pipiens_pangenome
# ============================================================

# --- 1. Install / load packages -------------------------------------------
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

# --- 2. Paths and constants -----------------------------------------------
proj <- "C:/Users/tam35/cx_pipiens_pangenome"   # change if running elsewhere
res_dir <- file.path(proj, "results")
out_dir <- file.path(proj, "figures")
anvio_dir <- file.path(proj, "anvio_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples_ingroup <- c("Cx_molestus", "Cx_pallens", "Cx_pipiens", "Cx_quinquefasciatus")
samples_all     <- c(samples_ingroup, "Cx_tarsalis")
species_pretty <- c(Cx_molestus = "Cx. molestus", Cx_pallens = "Cx. pallens",
                    Cx_pipiens = "Cx. pipiens", Cx_quinquefasciatus = "Cx. quinquefasciatus",
                    Cx_tarsalis = "Cx. tarsalis")

# ---- Unified colour palette (MDPI-friendly, color-blind safe) ----------------
# Each color category serves ONE semantic purpose throughout the paper so the
# reader builds a stable visual association. Distinct *families* avoid the
# trap of "orange = molestus AND CAFE-contraction".
#
#   Okabe-Ito (categorical, reserved for SPECIES only)
#   Sequential Blues  → pangenome core/shell/cloud (ordered: dense → sparse)
#   Diverging RdYlGn  → BUSCO completeness (ordered: best → worst)
#   Diverging G→R     → CAFE expansion / contraction (universal convention)
#   Brewer Set2       → TE composition categories (qualitative, distinct from
#                       species, distinct from compartment)
#   Brewer Set1 (3)   → chromosomes (qualitative, distinct from species and TE)
#
# Greys for neutral/missing/connector elements.

pal_grey   <- list(dark = "#444444", mid = "#999999", light = "#DDDDDD")

# Species — Okabe-Ito (color-blind safe categorical). Used ONLY for species.
species_colors <- c(
  Cx_molestus         = "#E69F00",   # orange
  Cx_pallens          = "#56B4E9",   # sky blue
  Cx_pipiens          = "#009E73",   # bluish green
  Cx_quinquefasciatus = "#CC79A7",   # reddish purple
  Cx_tarsalis         = "#999999"    # grey
)

# Pangenome compartment — sequential Blues (dense → sparse)
comp_colors <- c(
  core  = "#08519C",   # darkest
  shell = "#6BAED6",   # medium
  cloud = "#C6DBEF"    # lightest
)

# BUSCO completeness — RdYlGn diverging (best → worst)
busco_colors <- c(
  "Single copy" = "#1A9850",   # dark green (good)
  "Multi copy"  = "#A6D96A",   # light green
  "Fragmented"  = "#FDAE61",   # amber
  "Missing"     = "#D73027"    # red (bad)
)

# CAFE — green (expansion) / red (contraction); universal scientific convention
cafe_colors <- c(
  Expansion   = "#4DAC26",
  Contraction = "#CA0020"
)

# TE composition — Brewer Set2 (qualitative, distinct from species palette)
te_colors <- c(
  "Interspersed TE" = "#66C2A5",   # teal
  "Simple repeats"  = "#FC8D62",   # salmon
  "Low complexity"  = "#8DA0CB",   # lavender
  "Unmasked"        = "#B3B3B3"    # neutral grey
)

# Chromosomes — Brewer Set1 (qualitative, distinct again)
chrom_colors <- c(
  chr1 = "#377EB8",   # blue
  chr2 = "#4DAF4A",   # green
  chr3 = "#E41A1C"    # red
)

# Convenience aliases used in fig9 anvi'o UpSet (dots / connectors)
pal <- list(dark_grey = pal_grey$dark, light_grey = pal_grey$light)

# Shared theme: Arial sans-serif (MDPI Insects requires sans-serif), minimal grid
# MDPI specs: single-column figures ≤ 9 cm wide, double-column ≤ 18.5 cm,
# at least 300 dpi (we save at 300 dpi PNG + vector PDF).
theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(panel.grid.minor = element_blank(),
          axis.title  = element_text(size = base_size - 1),
          axis.text   = element_text(size = base_size - 2, color = "black"),
          plot.title  = element_text(size = base_size, face = "plain"),
          legend.title= element_text(size = base_size - 1, face = "bold"),
          legend.text = element_text(size = base_size - 1))
}

# Helper: save both PNG and PDF
save_fig <- function(p, name, width, height, dpi = 300) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height)
  message("  wrote ", name, ".{png,pdf}")
}

italic_sp <- function(s) bquote(italic(.(s)))

# =========================================================================
# Figure 1 — Assembly quality (BUSCO proteins + QUAST stats table)
# =========================================================================
fig1_assembly_quality <- function() {
  message("Figure 1: Assembly quality")
  # ---- BUSCO from JSON
  busco <- lapply(samples_all, function(s) {
    fp <- file.path(res_dir, "busco_proteins", s,
                    paste0("short_summary.specific.diptera_odb10.", s, ".json"))
    if (!file.exists(fp)) return(NULL)
    j <- jsonlite::fromJSON(fp)
    r <- j$results
    data.frame(sample = s,
               "Single copy" = r[["Single copy percentage"]],
               "Multi copy"  = r[["Multi copy percentage"]],
               "Fragmented"  = r[["Fragmented percentage"]],
               "Missing"     = r[["Missing percentage"]],
               check.names = FALSE)
  })
  busco_df <- do.call(rbind, busco)
  busco_long <- pivot_longer(busco_df, -sample,
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
         x = expression("BUSCO (% of 3,285 " * italic("diptera_odb10") * " markers)"), y = NULL,
         fill = NULL) +
    theme_pub() +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank())

  # ---- QUAST table (rendered as a tableGrob)
  quast <- fread(file.path(res_dir, "quast", "transposed_report.tsv"))
  quast$Assembly <- factor(quast$Assembly, levels = samples_all)
  tab <- quast[order(Assembly), .(
    Species = species_pretty[as.character(Assembly)],
    "Size (Mb)" = round(`Total length (>= 0 bp)` / 1e6, 1),
    "Contigs"  = `# contigs`,
    "N50 (Mb)" = round(N50 / 1e6, 1),
    "GC (%)"   = round(`GC (%)`, 2)
  )]
  tab_grob <- gridExtra::tableGrob(tab, rows = NULL,
    theme = gridExtra::ttheme_minimal(
      core    = list(fg_params = list(cex = 0.85, fontface = ifelse(grepl("^Cx", tab$Species), 3, 1))),
      colhead = list(fg_params = list(cex = 0.9, fontface = "bold"))))
  p1b <- wrap_elements(full = tab_grob) +
    labs(title = "(b) Assembly statistics") + theme_pub() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())

  out <- p1a + p1b + plot_layout(widths = c(1.3, 1))
  print(out)
  save_fig(out, "Figure_1_assembly_quality", width = 13, height = 5)
}

# =========================================================================
# Figure 2 — Pangenome composition (donut + per-form bars + form-specific cloud)
# =========================================================================
fig2_pangenome <- function() {
  message("Figure 2: Pangenome composition")
  summ <- fread(file.path(res_dir, "pangenome", "pangenome_summary.tsv"))
  parts <- fread(file.path(res_dir, "pangenome", "partitioned_orthogroups.tsv"))
  parts[, comp_main := sub("_Cx_.*_specific$", "", compartment)]

  # ---- (a) Donut with labels INSIDE the wedges (big slices) and outside for cloud
  main3 <- summ[compartment %in% c("core","shell","cloud")]
  main3 <- main3[match(c("core","shell","cloud"), compartment)]
  main3[, pct := n_orthogroups / sum(n_orthogroups) * 100]
  main3[, y_mid := cumsum(n_orthogroups) - n_orthogroups / 2]
  # Labels: big slices (>=10%) inside the wedge with white bold text;
  # tiny cloud slice (~5%) outside on a single line.
  main3[, label_inside := ifelse(pct >= 10,
                                  paste0(tools::toTitleCase(compartment), "\n",
                                         format(n_orthogroups, big.mark=","), "\n(",
                                         sprintf("%.1f", pct), "%)"),
                                  "")]
  main3[, label_outside := ifelse(pct < 10,
                                   paste0(tools::toTitleCase(compartment), " ",
                                          format(n_orthogroups, big.mark=","),
                                          " (", sprintf("%.1f", pct), "%)"),
                                   "")]
  main3$compartment <- factor(main3$compartment, levels = c("cloud","shell","core"))
  total_og <- sum(main3$n_orthogroups)

  p2a <- ggplot(main3, aes(x = 2, y = n_orthogroups, fill = compartment)) +
    geom_bar(stat = "identity", color = "white", width = 1) +
    coord_polar(theta = "y", start = 0) +
    xlim(0.4, 3.0) +
    scale_fill_manual(values = comp_colors, breaks = c("core","shell","cloud"),
                      labels = function(x) tools::toTitleCase(x)) +
    # Labels inside Core and Shell wedges (white text on dark/medium blue)
    geom_text(aes(x = 2.0, y = y_mid, label = label_inside),
              size = 3.1, lineheight = 0.95, color = "white",
              fontface = "bold", show.legend = FALSE) +
    # Label outside Cloud wedge (tiny slice, single line, dark text)
    geom_text(aes(x = 2.75, y = y_mid, label = label_outside),
              size = 3.0, color = "black", show.legend = FALSE) +
    annotate("text", x = 0.4, y = 0, label = format(total_og, big.mark=","),
             size = 6, fontface = "bold") +
    annotate("text", x = 0.4, y = 0, label = "orthogroups",
             size = 3.0, vjust = 3) +
    labs(title = "(a) Pangenome partition", fill = NULL) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(size = 11, hjust = 0.5, face = "plain"),
          legend.position = "none",
          plot.margin = margin(2, 2, 2, 2))

  # ---- (b) Per-form stacked gene counts
  per_form <- rbindlist(lapply(samples_ingroup, function(s) {
    data.table(species = s,
               core  = sum(parts[comp_main == "core",  ..s][[1]]),
               shell = sum(parts[comp_main == "shell", ..s][[1]]),
               cloud = sum(parts[comp_main == "cloud", ..s][[1]]))
  }))
  per_form_long <- melt(per_form, id.vars = "species",
                        variable.name = "compartment", value.name = "genes")
  per_form_long$compartment <- factor(per_form_long$compartment, levels = c("cloud","shell","core"))
  totals <- per_form_long[, .(total = sum(genes)), by = species]

  p2b <- ggplot(per_form_long, aes(x = species, y = genes, fill = compartment)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.4) +
    geom_text(data = totals, aes(x = species, y = total, label = format(total, big.mark=",")),
              vjust = -0.4, inherit.aes = FALSE, size = 3.2) +
    scale_fill_manual(values = comp_colors, breaks = c("core","shell","cloud"),
                      labels = function(x) tools::toTitleCase(x)) +
    scale_x_discrete(labels = function(x) parse(text = paste0("italic('", species_pretty[x], "')"))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = "(b) Genes per species by compartment", x = NULL, y = "Genes", fill = NULL) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          legend.position = "bottom")

  # ---- (c) Form-specific cloud
  cloud_specific <- summ[grepl("_specific$", compartment)]
  cloud_specific[, species := sub("cloud_(Cx_.*)_specific", "\\1", compartment)]
  cloud_specific$species <- factor(cloud_specific$species, levels = samples_ingroup)
  cloud_specific$species_color <- species_colors[as.character(cloud_specific$species)]

  p2c <- ggplot(cloud_specific, aes(x = species, y = n_orthogroups, fill = species)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.4, show.legend = FALSE) +
    geom_text(aes(label = n_orthogroups), vjust = -0.4, size = 3.5) +
    scale_fill_manual(values = species_colors) +
    scale_x_discrete(labels = function(x) parse(text = paste0("italic('", species_pretty[x], "')"))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "(c) Form-specific cloud orthogroups", x = NULL, y = "Orthogroups") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  out <- p2a + p2b + p2c + plot_layout(ncol = 3, widths = c(1.3, 1.1, 1.0))
  print(out)
  save_fig(out, "Figure_2_pangenome", width = 14.5, height = 5.2)
}

# =========================================================================
# Figure 3 — Species tree with concordance (ggtree)
# =========================================================================
fig3_species_tree <- function() {
  message("Figure 3: Species tree (ggtree)")
  # Read concord.cf.tree — bootstrap/gCF/sCF combined as a text node label
  tr <- ape::read.tree(file.path(res_dir, "phylo", "concord.cf.tree"))
  # ggtree wants species labels prettified
  tr$tip.label <- species_pretty[tr$tip.label]

  p <- ggtree(tr, size = 0.9) +
    geom_tiplab(fontface = "italic", offset = 0.002, size = 4.2) +
    geom_nodelab(aes(label = label), hjust = -0.10, vjust = -0.6,
                 fontface = "bold", size = 3.6, color = "#222222") +
    geom_treescale(x = 0, y = -0.4, width = 0.05, fontsize = 3.2,
                   linesize = 0.7, offset = 0.1) +
    xlim(0, max(ggtree(tr)$data$x) * 1.55) +
    labs(caption = paste0(
      "Node labels: ultrafast bootstrap / gene concordance / site concordance.\n",
      "Tree rooting is illustrative; Cx. tarsalis outgroup not in current run."
    )) +
    theme(plot.caption = element_text(hjust = 0.5, size = 9,
                                       face = "italic", color = "#444444"),
          plot.margin = margin(15, 30, 15, 15))

  print(p)
  save_fig(p, "Figure_3_species_tree", width = 9.5, height = 6)
}

# =========================================================================
# Figure 4 — CAFE gene family evolution
# =========================================================================
fig4_cafe <- function() {
  message("Figure 4: CAFE")
  clade <- fread(file.path(res_dir, "cafe", "output", "Gamma_clade_results.txt"))
  # Strip <N> branch tags from Taxon_ID
  clade[, species := sub("<[0-9]+>$", "", `#Taxon_ID`)]
  clade$species <- factor(clade$species, levels = samples_ingroup)

  long <- melt(clade[, .(species, Expansion = Increase, Contraction = Decrease)],
               id.vars = "species", variable.name = "type", value.name = "n")
  long$type <- factor(long$type, levels = c("Expansion","Contraction"))

  # Pull lambda, alpha, n_sig for subtitle
  res_lines <- readLines(file.path(res_dir, "cafe", "output", "Gamma_results.txt"))
  lambda <- as.numeric(sub("Lambda: ", "", grep("^Lambda:", res_lines, value = TRUE)))
  alpha  <- as.numeric(sub("Alpha: ",  "", grep("^Alpha:",  res_lines, value = TRUE)))
  n_sig <- length(readLines(file.path(res_dir, "cafe", "significant_families.tsv"))) - 1
  # Use plotmath expression — Windows R can't render lambda/alpha/superscripts via Unicode
  subtitle_expr <- bquote(.(n_sig) ~ "significant families (p<0.05);  " *
                          lambda == .(sprintf("%.4f", lambda)) ~ "," ~
                          alpha == .(sprintf("%.3f", alpha)) *
                          ";  binomial p<2x10^-16 each branch")

  p <- ggplot(long, aes(x = species, y = n, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.78),
             width = 0.7, color = "white", linewidth = 0.4) +
    geom_text(aes(label = format(n, big.mark = ",")),
              position = position_dodge(width = 0.78), vjust = -0.4, size = 3.4) +
    scale_fill_manual(values = cafe_colors, name = NULL) +
    scale_x_discrete(labels = function(x) parse(text = paste0("italic('", species_pretty[x], "')"))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "Gene families", subtitle = subtitle_expr) +
    theme_pub() +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 9, color = "#444444"))

  print(p)
  save_fig(p, "Figure_4_cafe", width = 10, height = 5.5)
}

# =========================================================================
# Figure 5 — Key gene families heatmap
# =========================================================================
fig5_key_families <- function() {
  message("Figure 5: Key gene families")
  d <- fread(file.path(res_dir, "functional", "key_families_wide.tsv"))
  # Order families by mean copy number (highest first); rowMeans needs a matrix
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
    scale_y_discrete(limits = rev) +   # highest-mean family (P450) on top
    labs(x = NULL, y = "Gene family") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          panel.grid = element_blank())

  print(p)
  save_fig(p, "Figure_9_key_families", width = 8.5, height = 5.5)
}

# =========================================================================
# Figure 6 — TE composition
# =========================================================================
fig6_te_composition <- function() {
  message("Figure 6: TE composition")
  rows <- lapply(samples_ingroup, function(s) {
    tbl_lines <- readLines(file.path(res_dir, "repeats", s, paste0(s, ".fasta.tbl")))
    total      <- as.numeric(sub(".*total length:\\s+(\\d+)\\s+bp.*", "\\1",
                                  grep("^total length:", tbl_lines, value = TRUE)[1]))
    interspersed <- as.numeric(sub(".*Total interspersed repeats:\\s+(\\d+)\\s+bp.*", "\\1",
                                    grep("Total interspersed repeats:", tbl_lines, value = TRUE)[1]))
    simple     <- as.numeric(sub(".*Simple repeats:\\s+\\d+\\s+(\\d+)\\s+bp.*", "\\1",
                                  grep("^Simple repeats:", tbl_lines, value = TRUE)[1]))
    low        <- as.numeric(sub(".*Low complexity:\\s+\\d+\\s+(\\d+)\\s+bp.*", "\\1",
                                  grep("^Low complexity:", tbl_lines, value = TRUE)[1]))
    unmasked <- total - interspersed - simple - low
    data.table(sample = s,
               "Interspersed TE" = interspersed/total*100,
               "Simple repeats"  = simple/total*100,
               "Low complexity"  = low/total*100,
               "Unmasked"        = unmasked/total*100,
               check.names = FALSE)
  })
  d <- rbindlist(rows)
  long <- melt(d, id.vars = "sample", variable.name = "category", value.name = "pct")
  # For horizontal stacked bars ggplot stacks factor levels in reverse, so to
  # put Interspersed TE on the LEFT we keep natural factor order and use
  # position_stack(reverse = TRUE) which flips the rendering once.
  long$category <- factor(long$category,
                          levels = c("Interspersed TE","Simple repeats","Low complexity","Unmasked"))
  long$sample <- factor(long$sample, levels = rev(samples_ingroup))

  p <- ggplot(long, aes(x = pct, y = sample, fill = category)) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE),
             color = "white", linewidth = 0.3) +
    geom_text(data = long[long$pct > 4, ],
              aes(label = sprintf("%.1f%%", pct)),
              position = position_stack(vjust = 0.5, reverse = TRUE), size = 3.2,
              color = ifelse(long$category[long$pct > 4] == "Interspersed TE",
                             "white", "black")) +
    scale_fill_manual(values = te_colors,
                      breaks = c("Interspersed TE","Simple repeats","Low complexity","Unmasked")) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 100.5)) +
    scale_y_discrete(labels = function(x) parse(text = paste0("italic('", species_pretty[x], "')"))) +
    labs(x = "Genome fraction (%)", y = NULL, fill = NULL) +
    theme_pub() +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank())

  print(p)
  save_fig(p, "Figure_7_te_composition", width = 10, height = 4.8)
}

# =========================================================================
# Figure 7 — Pairwise whole-genome ANI heatmap
# =========================================================================
fig7_ani <- function() {
  message("Figure 7: ANI heatmap")
  df <- fread(file.path(res_dir, "synteny", "ani_matrix.tsv"))
  setnames(df, "sample", "row_sample")
  # Restrict to ingroup (tarsalis is NA)
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
    coord_equal() +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          panel.grid = element_blank())

  print(p)
  save_fig(p, "Figure_5_ani", width = 7, height = 5.5)
}

# =========================================================================
# Figure 8 — Pairwise synteny dotplots
# =========================================================================
fig8_synteny <- function() {
  message("Figure 8: Synteny")
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
    n_skip <- count_header_lines(fp)
    d <- fread(fp, sep = "\t", header = FALSE, fill = TRUE, skip = n_skip,
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

    ggplot(d, aes(x = x_pos/1e6, y = y_pos/1e6, color = x_label)) +
      geom_point(size = 0.35, alpha = 0.55) +
      scale_color_manual(values = chrom_colors,
                         name = "Reference chromosome",
                         guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(title = paste0("(", panel_letter, ") ", species_pretty[sx], " vs ", species_pretty[sy]),
           x = paste0(species_pretty[sx], " (Mb)"),
           y = paste0(species_pretty[sy], " (Mb)")) +
      theme_pub(base_size = 10) +
      theme(axis.title = element_text(face = "italic", size = 9))
  }

  pairs_list <- list(
    c("Cx_pallens", "Cx_pipiens"),
    c("Cx_pallens", "Cx_quinquefasciatus"),
    c("Cx_molestus", "Cx_pallens"),
    c("Cx_molestus", "Cx_pipiens"),
    c("Cx_pipiens", "Cx_quinquefasciatus"),
    c("Cx_molestus", "Cx_quinquefasciatus"))
  plots <- mapply(function(p, l) plot_pair(p[1], p[2], l),
                  pairs_list, letters[seq_along(pairs_list)], SIMPLIFY = FALSE)
  out <- wrap_plots(plots, ncol = 3, nrow = 2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  print(out)
  save_fig(out, "Figure_6_synteny", width = 13, height = 8.5)
}

# =========================================================================
# Figure 9 — Anvi'o pangenome composition (UpSet + per-genome + size dist)
# =========================================================================
fig9_anvio <- function() {
  message("Figure 9: Anvi'o pangenome")
  fp <- file.path(anvio_dir, "Culex_pipiens_complex_gene_clusters_summary.txt")
  d <- fread(fp, select = c("gene_cluster_id", "genome_name"))
  pa <- dcast(d[, .N, by = .(gene_cluster_id, genome_name)],
              gene_cluster_id ~ genome_name, value.var = "N", fill = 0)
  for (s in samples_ingroup) if (!s %in% names(pa)) pa[[s]] <- 0
  pa <- pa[, c("gene_cluster_id", samples_ingroup), with = FALSE]
  for (s in samples_ingroup) pa[[s]] <- as.integer(pa[[s]] > 0)
  pa[, pattern := apply(.SD, 1, paste, collapse = ""), .SDcols = samples_ingroup]
  pat_counts <- pa[, .N, by = pattern][order(-N)]
  pat_counts[, x := .I]
  total_clusters <- nrow(pa)

  # (a) UpSet-style intersection bars
  p_bars <- ggplot(pat_counts, aes(x = x, y = N)) +
    geom_bar(stat = "identity", fill = pal$dark_grey, color = "white", linewidth = 0.3) +
    geom_text(aes(label = format(N, big.mark = ",")), vjust = -0.4, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "Gene clusters", title = "(a) Cluster intersections") +
    theme_pub() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank())

  # Matrix below bars: dots showing membership
  mat_rows <- rbindlist(lapply(seq_len(nrow(pat_counts)), function(i) {
    bits <- as.integer(strsplit(pat_counts$pattern[i], "")[[1]])
    data.table(x = i,
               species = samples_ingroup,
               present = bits)
  }))
  mat_rows$species <- factor(mat_rows$species, levels = rev(samples_ingroup))
  mat_rows$y <- as.integer(mat_rows$species)
  # Vertical connectors per intersection (lowest to highest "present" y)
  conns <- mat_rows[present == 1, .(ymin = min(y), ymax = max(y)), by = x][ymin < ymax]

  p_mat <- ggplot() +
    geom_segment(data = conns, aes(x = x, xend = x, y = ymin, yend = ymax),
                 color = pal$dark_grey, linewidth = 0.8) +
    geom_point(data = mat_rows[present == 1, ],
               aes(x = x, y = y, color = species), size = 3.6) +
    geom_point(data = mat_rows[present == 0, ],
               aes(x = x, y = y), color = pal$light_grey, size = 3.6) +
    scale_color_manual(values = species_colors, guide = "none") +
    scale_y_continuous(breaks = seq_along(samples_ingroup),
                       labels = function(b) parse(text = paste0("italic('", species_pretty[rev(samples_ingroup)[b]], "')"))) +
    scale_x_continuous(breaks = NULL) +
    labs(x = "Intersection (sorted by size, left = largest)", y = NULL) +
    theme_void(base_size = 10) +
    theme(axis.text.y = element_text(hjust = 1),
          axis.title.x = element_text(size = 10),
          plot.margin = margin(2, 5, 5, 5))

  # (b) Per-genome cluster count
  per_g <- data.table(species = factor(samples_ingroup, levels = samples_ingroup),
                      n = sapply(samples_ingroup, function(s) sum(pa[[s]])))

  p_pg <- ggplot(per_g, aes(x = n, y = species, fill = species)) +
  p_pg <- ggplot(per_g, aes(x = n, y = species, fill = species)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.3, show.legend = FALSE) +
    geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.15, size = 3.3) +
    scale_fill_manual(values = species_colors) +
    scale_y_discrete(labels = function(x) parse(text = paste0("italic('", species_pretty[x], "')")),
                     limits = rev) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = "Clusters per genome", y = NULL, title = "(b) Per-genome cluster count") +
    theme_pub()

  # (c) Cluster size distribution
  sizes <- d[, .N, by = gene_cluster_id]
  p_sz <- ggplot(sizes, aes(x = N)) +
    geom_histogram(bins = 30, fill = pal$dark_grey, color = "white", linewidth = 0.3) +
    scale_x_log10() +
    annotation_logticks(sides = "b", size = 0.4) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("median %d | max %d", median(sizes$N), max(sizes$N)),
             hjust = 1.05, vjust = 1.5, size = 3.3, color = "#444") +
    labs(x = "Genes per cluster", y = "# clusters", title = "(c) Cluster-size distribution") +
    theme_pub()

  # Layout: left = (a)+matrix; right = (b)/(c)
  left_col <- p_bars / p_mat + plot_layout(heights = c(3.5, 1.5))
  right_col <- p_pg / p_sz + plot_layout(heights = c(1, 1))
  out <- left_col | right_col
  out <- out + plot_layout(widths = c(2.3, 1))
  print(out)
  save_fig(out, "Figure_8_anvio_pangenome", width = 14, height = 8.5)
}

# =========================================================================
# Main
# =========================================================================
main <- function() {
  cat("==> Generating all figures into:", out_dir, "\n\n")
  fig1_assembly_quality()
  fig2_pangenome()
  fig3_species_tree()
  fig4_cafe()
  fig5_key_families()
  fig6_te_composition()
  fig7_ani()
  fig8_synteny()
  fig9_anvio()
  cat("\n==> All 9 figures written.\n")
}

# Run only if executed directly (not on source())
if (sys.nframe() == 0) main()
