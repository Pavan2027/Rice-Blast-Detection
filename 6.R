# =====================================================================
# SCRIPT 6: Heatmap Visualizations
# Rice Blast Resistance Project — VIT BCSE207L
# Continuation of Scripts 1–5
# Generates:
#   heatmap_1_grm.png        — Genomic Relationship Matrix
#   heatmap_2_snp_genotype.png — Top GWAS SNPs × Plants genotype pattern
# =====================================================================

load("rice_blast_gwas.RData")    # gives: X_encoded, Y_ml, gwas_results, top_snps, pca_res
load("rice_blast_final.RData")   # gives: gblup_fit, lgb_model, tabnet_model, summary_table

# Install pheatmap if not already installed
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos = "https://cloud.r-project.org")
}

library(pheatmap)
library(ggplot2)
library(rrBLUP)

# =====================================================================
# HEATMAP 1: Genomic Relationship Matrix (GRM)
# Shows genetic similarity/kinship between all 364 rice accessions.
# GBLUP uses this matrix internally — visualizing it explains HOW
# GBLUP works and why it captures population structure.
# =====================================================================
cat("\n[1/2] Building Genomic Relationship Matrix for heatmap...\n")
cat("      (This may take ~30 seconds)\n\n")

# Build G matrix from all 35,800 SNPs — same as Script 4
G_matrix <- A.mat(X_encoded)

# Sort accessions by blast score so resistant plants cluster together
blast_order <- order(Y_ml$Blast_Score)
G_sorted    <- G_matrix[blast_order, blast_order]

# Annotation bar: disease status for each accession (left side of heatmap)
status_sorted <- Y_ml$Disease_Status[blast_order]
row_annotation <- data.frame(
  Disease_Status = status_sorted,
  row.names      = rownames(G_sorted)
)

# Annotation bar: subpopulation
subpop_sorted <- Y_ml$Subpopulations[blast_order]
subpop_sorted[subpop_sorted == "Not available"] <- "Unknown"
row_annotation$Subpopulation <- subpop_sorted

# Build Subpopulation colors dynamically from what's actually in the data
# This prevents the "factor levels do not match" error regardless of dataset
subpop_palette <- c(
  "IND"      = "#3498db",
  "TEJ"      = "#9b59b6",
  "TRJ"      = "#e67e22",
  "AUS"      = "#1abc9c",
  "AROMATIC" = "#f39c12",
  "ADMIX"    = "#e91e63",
  "Unknown"  = "#bdc3c7"
)

# Get the subpopulations that actually appear in THIS dataset
actual_subpops <- unique(subpop_sorted)

# Only keep colors for subpopulations that exist — pheatmap requires exact match
subpop_colors_used <- subpop_palette[names(subpop_palette) %in% actual_subpops]

# If any subpopulation in data has no assigned color, assign gray as fallback
missing_subpops <- setdiff(actual_subpops, names(subpop_palette))
if (length(missing_subpops) > 0) {
  fallback <- setNames(rep("#95a5a6", length(missing_subpops)), missing_subpops)
  subpop_colors_used <- c(subpop_colors_used, fallback)
}

annotation_colors <- list(
  Disease_Status = c(
    "Resistant"   = "#27ae60",
    "Susceptible" = "#e74c3c"
  ),
  Subpopulation = subpop_colors_used
)

# Color palette for GRM values:
# Blue = genetically similar (high kinship), White = average, Red = dissimilar
grm_colors <- colorRampPalette(c("#2980b9", "#ecf0f1", "#c0392b"))(100)

cat("Saving GRM heatmap...\n")

pheatmap(
  G_sorted,
  color             = grm_colors,
  annotation_row    = row_annotation,
  annotation_col    = row_annotation,
  annotation_colors = annotation_colors,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  border_color      = NA,
  main              = "Genomic Relationship Matrix (GRM) — 364 Rice Accessions\nSorted by Blast Resistance Score | Blue = Genetically Similar",
  fontsize          = 11,
  legend_breaks     = c(-0.2, 0, 0.5, 1.0, 1.5),
  legend_labels     = c("-0.2", "0 (unrelated)", "0.5", "1.0 (self)", "1.5"),
  annotation_legend = TRUE,
  filename          = "heatmap_1_grm.png",
  width             = 14,
  height            = 12
)

file_size <- file.info("heatmap_1_grm.png")$size
cat("  -> heatmap_1_grm.png size:", file_size, "bytes\n")
cat("  -> heatmap_1_grm.png saved!\n\n")
# =====================================================================
# HEATMAP 2: Top GWAS SNPs × Plants Genotype Pattern
# =====================================================================
cat("[2/2] Building SNP genotype heatmap for top GWAS SNPs...\n\n")

# Extract the 19 GWAS-selected SNPs for all 364 plants
X_top <- X_encoded[, top_snps]

# Sort plants by blast score (most resistant at top/left)
blast_order_2  <- order(Y_ml$Blast_Score)
X_top_sorted   <- X_top[blast_order_2, ]
status_sorted_2 <- as.character(Y_ml$Disease_Status[blast_order_2])
blast_sorted    <- Y_ml$Blast_Score[blast_order_2]

# ---- Plant annotation (columns) ----
# Use simple string labels with NO special characters to avoid factor mismatch
blast_tier <- ifelse(blast_sorted <= 2, "Low 0-2",
                     ifelse(blast_sorted <= 5, "Medium 3-5", "High 6-9"))

plant_ann <- data.frame(
  Disease_Status = status_sorted_2,
  Blast_Tier     = blast_tier,
  row.names      = rownames(X_top_sorted),
  stringsAsFactors = FALSE
)

# ---- SNP annotation (rows) ----
# Match significance to SNPs — use simple labels only
gwas_top_sub <- gwas_results[gwas_results$SNP %in% top_snps, ]

sig_label <- ifelse(gwas_top_sub$neglog10P >= 5, "Very Strong p<1e-5",
                    ifelse(gwas_top_sub$neglog10P >= 4, "Strong p<1e-4",
                           "Suggestive"))

snp_ann <- data.frame(
  Significance = sig_label[match(colnames(X_top_sorted), gwas_top_sub$SNP)],
  row.names    = colnames(X_top_sorted),
  stringsAsFactors = FALSE
)

# Replace any NA significance with "Suggestive" as fallback
snp_ann$Significance[is.na(snp_ann$Significance)] <- "Suggestive"

# ---- Verify all factor levels match annotation_colors BEFORE plotting ----
cat("Disease_Status levels in data:", unique(plant_ann$Disease_Status), "\n")
cat("Blast_Tier levels in data:",     unique(plant_ann$Blast_Tier), "\n")
cat("Significance levels in data:",   unique(snp_ann$Significance), "\n")

# ---- Annotation colors — keys must EXACTLY match labels above ----
snp_annotation_colors <- list(
  Disease_Status = c(
    "Resistant"   = "#27ae60",
    "Susceptible" = "#e74c3c"
  ),
  Blast_Tier = c(
    "Low 0-2"    = "#1a7a4a",
    "Medium 3-5" = "#f39c12",
    "High 6-9"   = "#c0392b"
  ),
  Significance = c(
    "Suggestive"        = "#f1c40f",
    "Strong p<1e-4"     = "#e67e22",
    "Very Strong p<1e-5" = "#c0392b"
  )
)

# ---- Clean SNP labels ----
snp_labels <- gsub("^id|^ud|^wd", "", colnames(X_top_sorted))

# ---- Gap line position ----
n_resistant <- sum(status_sorted_2 == "Resistant")
cat("Gap will be placed after column:", n_resistant, "\n")

# ---- Plot ----
cat("Saving SNP genotype heatmap...\n")

pheatmap(
  t(X_top_sorted),
  color             = snp_colors,
  annotation_col    = plant_ann,
  annotation_row    = snp_ann,
  annotation_colors = snp_annotation_colors,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  labels_row        = snp_labels,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  border_color      = NA,
  main = paste0(
    "Genotype Pattern of Top 19 GWAS-Significant SNPs\n",
    "Plants sorted by Blast Score ",
    "(Left = Resistant -> Right = Susceptible)\n",
    "Yellow = Reference Allele | Blue = Alternate Allele"
  ),
  fontsize          = 11,
  fontsize_row      = 9,
  legend_breaks     = c(0, 2),
  legend_labels     = c("Reference (0)", "Alternate (2)"),
  gaps_col          = n_resistant,
  filename          = "heatmap_2_snp_genotype.png",  # pheatmap saves directly
  width             = 10,   # inches
  height            = 14    # inches
)

# Verify
file_size <- file.info("heatmap_2_snp_genotype.png")$size
cat("  -> heatmap_2_snp_genotype.png size:", file_size, "bytes\n")
if (file_size < 10000) {
  cat("  WARNING: Still small — something went wrong\n")
} else {
  cat("  -> SNP genotype heatmap saved successfully!\n\n")
}
dev.off()

# Verify file was written properly
file_size <- file.info("heatmap_2_snp_genotype.png")$size
cat("  -> heatmap_2_snp_genotype.png size:", file_size, "bytes\n")
if (file_size < 10000) {
  cat("  WARNING: File is suspiciously small — check error above\n")
} else {
  cat("  -> SNP genotype heatmap saved successfully!\n\n")
}

# =====================================================================
# BONUS: Print interpretation summary to console
# =====================================================================
cat("==========================================\n")
cat("          HEATMAP INTERPRETATION          \n")
cat("==========================================\n\n")

cat("Heatmap 1 — GRM:\n")
cat("  Diagonal = 1.0 (each plant compared to itself)\n")
cat("  Off-diagonal clusters = genetically related subpopulations\n")
cat("  Resistant plants (green) appear as distinct genomic clusters\n")
cat("  This structure is what GBLUP models, explaining its high sensitivity\n\n")

cat("Heatmap 2 — SNP Genotype Pattern:\n")
cat("  Each row = one of the 19 GWAS-significant SNPs\n")
cat("  Each column = one of 364 rice accessions\n")
cat("  Yellow = carries reference allele | Blue = carries alternate allele\n")
cat("  Vertical gap separates Resistant (left) from Susceptible (right)\n")
cat("  SNPs with clear color shifts at the gap = strongest blast predictors\n\n")

# Print top 5 most significant SNPs for reference
cat("Top 5 GWAS SNPs (for report):\n")
top5 <- head(gwas_results[order(gwas_results$Pval), ], 5)
top5$Pval_formatted <- formatC(top5$Pval, format = "e", digits = 2)
print(top5[, c("SNP", "Beta", "Pval_formatted", "neglog10P")], row.names = FALSE)

cat("\n==========================================\n")
cat("DONE! Files saved:\n")
cat("  heatmap_1_grm.png\n")
cat("  heatmap_2_snp_genotype.png\n")
cat("==========================================\n")
