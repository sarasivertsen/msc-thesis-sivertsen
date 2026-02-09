# README
# How to run:
#   Rscript: count_table_analysis.R, dataset: fd_count.csv
# Workflow: 
#   1. CA biplot 
#   2. Hierarchical clustering analysis 
#   3. Fisher's exact test 
# Notes:
#   - The script installs necessary packages if needed 
#   - Outputs are written to an output folder

# ----- Installing packages, setting up directories, defining datasets, storing results ----- #

# Installing packages
required_pkgs <- c("FactoMineR", "ggplot2", "ggrepel", "scales")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}

# Creating directories and output folder
data_dir <- "flavor_data"
results_dir <- "flavor_results"

if(!dir.exists(data_dir)) {
  stop("Data directory does not exist: ", normalizePath(data_dir, mustWork = FALSE))
}

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

out_path <- function(name) file.path(results_dir, name)

# Loading input data
count_file <- file.path(data_dir, "fd_count.csv")
if (!file.exists(count_file)) {
  stop("Data file not found ", count_file)
}

count_df <- read.csv(count_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!"additive" %in% names(count_df)) {
  stop("Expected column 'additive' not in fd_count.csv")
}

# Building the count matrix 
flavor_names <- setdiff(names(count_df), "additive")
counts_mat <- as.matrix(count_df[, flavor_names])
rownames(counts_mat) <- count_df$additive

if (!all(sapply(counts_mat, is.numeric))) stop("All flavor columns must be numeric counts.")

# ---------- Correspondence Analysis ---------- #
facto_res <- FactoMineR::CA(counts_mat, graph = FALSE)

ca_rowcoord <- data.frame(
  additive = rownames(counts_mat), 
  facto_res$row$coord, 
  row.names = NULL
)

ca_colcoord <- data.frame(
  flavor = colnames(counts_mat), 
  facto_res$col$coord, 
  row.names = NULL
)

ca_plot <- ggplot2::ggplot() +
  ggplot2::geom_point(data = ca_rowcoord,
                      ggplot2::aes(x = Dim.1, y = Dim.2),
                      color = "navy", size = 2, shape = 19) +
  ggrepel::geom_text_repel(data = ca_rowcoord,
                           ggplot2::aes(x = Dim.1, y = Dim.2, label = additive),
                           size = 4,
                           fontface = "bold",
                           color = "navy",
                           box.padding = 0.8,
                           point.padding = 0.2,
                           segment.color = NA,  
                           max.overlaps = Inf) +
  ggplot2::geom_point(data = ca_colcoord,
                      ggplot2::aes(x = Dim.1, y = Dim.2),
                      color = "firebrick", size = 2, shape = 17)  +
  ggrepel::geom_text_repel(data = ca_colcoord,
                           ggplot2::aes(x = Dim.1, y = Dim.2, label = flavor),
                           size = 3.5,
                           color = "firebrick",
                           segment.color = "grey70",
                           segment.size = 0.5,
                           box.padding = 0.3,
                           point.padding = 0.4,
                           force = 6,
                           max.overlaps = Inf) +
  ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = "dashed", linewidth = 0.2) +
  ggplot2::geom_vline(xintercept = 0, color = "grey50", linetype = "dashed", linewidth = 0.2) +
  ggplot2::labs(x = sprintf("CA Dim 1 (%.1f%%)", facto_res$eig[1, "percentage of variance"]),
                y = sprintf("CA Dim 2 (%.1f%%)", facto_res$eig[2, "percentage of variance"])) +
  ggplot2::theme_minimal() + 
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(), 
    panel.border = ggplot2::element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  )

# Save the biplot
png(out_path("ca_biplot.png"), width = 1200, height = 900, res = 150)
print(ca_plot)
dev.off()

# ----- Hierarchical clustering analysis (Bray-Curtis distances, avg. linkage) ----- #
bray_additives <- vegan::vegdist(counts_mat, method = "bray")
hc_additives <- hclust(bray_additives, method = "average")

png(out_path("hc_additives_dendrogram.png"), width = 900, height = 600, res = 150)
op <- par(no.readonly = TRUE)
par(mar = c(5, 5, 4, 2) + 0.1)  

plot(hc_additives,
     labels = rownames(counts_mat),
     main = "Hierarchical clustering: Flavor compounds",
     sub = "", xlab = "", ylab = "Bray-Curtis dissimilarity",
     hang = -0.1,
     ylim = c(0, 1),
     yaxt = "n")
axis(2, at = seq(0, 1, by = 0.1), las = 1)

par(op)
dev.off()

# ----- Fisher's exact tests (pairwise flavor x additive) ----- #
total_per_additive <- rowSums(counts_mat)
pair_tests <- list()
for (idx in seq_along(flavor_names)) {
  flav <- flavor_names[idx]
  for (a in seq_len(nrow(counts_mat))) {
    add_name <- rownames(counts_mat)[a]
    flav_in_add <- counts_mat[a, idx]
    flav_in_other_adds <- sum(counts_mat[-a, idx])
    notflav_in_add <- total_per_additive[a] - flav_in_add
    notflav_in_other <- sum(total_per_additive[-a]) - flav_in_other_adds
    tbl2 <- matrix(c(flav_in_add, flav_in_other_adds,
                     notflav_in_add, notflav_in_other),
                   nrow = 2, byrow = TRUE)
    ft2 <- fisher.test(tbl2)
    odds <- ft2$estimate
    if (is.null(odds) || length(odds) == 0) odds <- NA_real_
    pair_tests[[length(pair_tests) + 1]] <- data.frame(
      flavor = flav,
      additive = add_name,
      count_in_additive = flav_in_add,
      count_other_additives = flav_in_other_adds,
      p_value = ft2$p.value,
      odds_ratio = unname(odds),
      stringsAsFactors = FALSE
    )
  }
}

pair_df <- do.call(rbind, pair_tests)
pair_df$p_adj_BH <- p.adjust(pair_df$p_value, method = "BH")
pair_df <- pair_df[order(pair_df$p_adj_BH), ]

# Writing output files 
write.csv(ca_rowcoord, out_path("ca_rowcoord.csv"), row.names = FALSE)
write.csv(ca_colcoord, out_path("ca_colcoord.csv"), row.names = FALSE)
write.csv(pair_df, out_path("fisher_flavor_additive_pairwise.csv"), row.names = FALSE)
saveRDS(hc_additives, file = out_path("hclust_additives.rds"))

cat("\nOutputs written:\n",
    "- ", out_path("ca_rowcoord.csv"), ", ", out_path("ca_colcoord.csv"), " (CA scores)\n",
    "- ", out_path("ca_biplot.png"), " (CA biplot of additives and flavors)\n",
    "- ", out_path("hc_additives_dendrogram.png"), "(Bray-Curtis distance, average linkage)\n",
    "- ", out_path("fisher_flavor_additive_pairwise.csv"), " (per-flavor, per-additive 2x2 Fisher with BH)\n",
    "- ", out_path("hclust_additives.rds"), " (hclust objects)\n",
    sep = "")

# ----- Session information ----- #
writeLines(
  capture.output(sessionInfo()),
  out_path("sessionInfo_count.txt")
)