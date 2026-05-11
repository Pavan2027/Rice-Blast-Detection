# =====================================================================
# SCRIPT 3: PCA + GWAS
# =====================================================================

load("rice_blast_encoded.RData")
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is not installed. Install it once with: install.packages('ggplot2', repos='https://cloud.r-project.org')")
}
library(ggplot2)

# Load SNP annotation used for Manhattan plot positions
snp_info_file <- "data/RiceDiversity.44K.MSU6.SNP_Information.txt"
if (!file.exists(snp_info_file)) {
  stop(paste("SNP annotation file not found:", snp_info_file))
}
snp_info <- read.table(
  snp_info_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  fill = TRUE
)
# Drop trailing blank column if present in source file.
snp_info <- snp_info[, colnames(snp_info) != "", drop = FALSE]

required_snp_cols <- c("SNPID", "CHR", "corrected.MSU.v6")
missing_snp_cols <- setdiff(required_snp_cols, colnames(snp_info))
if (length(missing_snp_cols) > 0) {
  stop(paste("Missing required columns in SNP annotation:", paste(missing_snp_cols, collapse = ", ")))
}

# =====================================================================
# PART A: PCA — Population Structure
# =====================================================================
print("Running PCA... may take 30 seconds")

pca_res <- prcomp(X_encoded, center = TRUE, scale. = FALSE)

# How much variance does each PC explain?
variance_explained <- round(summary(pca_res)$importance[2, 1:10] * 100, 2)
print("Variance explained by top 10 PCs:")
print(variance_explained)

# Build PCA dataframe for plotting
pca_df <- as.data.frame(pca_res$x[, 1:5])
pca_df$SampleID <- rownames(pca_df)
pca_df$Disease_Status <- Y_ml$Disease_Status

# Plot PC1 vs PC2 colored by disease status
ggplot(pca_df, aes(x = PC1, y = PC2, color = Disease_Status)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("Resistant" = "#27ae60", "Susceptible" = "#e74c3c")) +
  labs(title = "PCA of 35,800 SNPs",
       subtitle = "Colored by Blast Disease Status",
       x = paste0("PC1 (", variance_explained[1], "%)"),
       y = paste0("PC2 (", variance_explained[2], "%)")) +
  theme_minimal()

# Save the plot
ggsave("pca_plot.png", width = 8, height = 6, dpi = 150)
print("PCA plot saved as pca_plot.png")

# =====================================================================
# PART B: GWAS — Find SNPs Associated with Blast Resistance
# =====================================================================
print("Running GWAS... this will take 3-5 minutes for 35,800 SNPs")

# Use PC1, PC2, PC3 as covariates to control for population structure
pca_covs <- pca_res$x[, 1:3]

# Binary outcome: 1 = Susceptible, 0 = Resistant
blast_binary <- as.integer(Y_ml$Disease_Status == "Susceptible")

# Run logistic regression for each SNP
n_snps <- ncol(X_encoded)
gwas_results <- data.frame(
  SNP  = colnames(X_encoded),
  Beta = NA_real_,
  Pval = NA_real_
)

for (i in seq_len(n_snps)) {
  snp <- X_encoded[, i]
  tryCatch({
    fit <- glm(blast_binary ~ snp + pca_covs[,1] + pca_covs[,2] + pca_covs[,3],
               family = binomial())
    coef_s <- summary(fit)$coefficients
    if ("snp" %in% rownames(coef_s)) {
      gwas_results$Beta[i] <- coef_s["snp", "Estimate"]
      gwas_results$Pval[i] <- coef_s["snp", "Pr(>|z|)"]
    }
  }, error = function(e) NULL)
  
  # Progress update every 5000 SNPs
  if (i %% 5000 == 0) print(paste("Processed", i, "of", n_snps, "SNPs..."))
}

# Remove SNPs that failed
gwas_results <- gwas_results[!is.na(gwas_results$Pval), ]
gwas_results$neglog10P <- -log10(gwas_results$Pval)

print(paste("GWAS complete.", nrow(gwas_results), "SNPs tested"))
print("Top 10 most significant SNPs:")
print(head(gwas_results[order(gwas_results$Pval), ], 10))
# Fix the merge — use SNPID column explicitly
gwas_plot <- merge(gwas_results, snp_info, by.x = "SNP", by.y = "SNPID")
if (nrow(gwas_plot) == 0) {
  stop("No SNPs matched between GWAS results and SNP annotation file.")
}
print(paste("SNPs with position info:", nrow(gwas_plot)))
# =====================================================================
# Manhattan Plot
# =====================================================================
gwas_plot$CHR <- suppressWarnings(as.integer(gwas_plot$CHR))
gwas_plot$corrected.MSU.v6 <- suppressWarnings(as.numeric(gwas_plot$corrected.MSU.v6))
gwas_plot <- gwas_plot[!is.na(gwas_plot$CHR) & !is.na(gwas_plot$corrected.MSU.v6), ]
gwas_plot <- gwas_plot[order(gwas_plot$CHR, gwas_plot$corrected.MSU.v6), ]

# Cumulative chromosome positions
chrom_offsets <- tapply(gwas_plot$corrected.MSU.v6, gwas_plot$CHR, max)
chrom_offsets <- cumsum(c(0, chrom_offsets[-length(chrom_offsets)]))
names(chrom_offsets) <- names(table(gwas_plot$CHR))

gwas_plot$bp_cumulative <- gwas_plot$corrected.MSU.v6 +
                            chrom_offsets[as.character(gwas_plot$CHR)]

# X-axis chromosome labels
axis_df <- aggregate(bp_cumulative ~ CHR, data = gwas_plot, FUN = mean)

# Plot
manhattan_plot <- ggplot(gwas_plot, aes(x = bp_cumulative, y = neglog10P,
                                         color = as.factor(CHR %% 2))) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("0" = "#2980b9", "1" = "#27ae60"), guide = "none") +
  geom_hline(yintercept = -log10(1e-4), linetype = "dashed", color = "orange", linewidth = 0.7) +
  geom_hline(yintercept = -log10(1e-6), linetype = "dashed", color = "red", linewidth = 0.7) +
  scale_x_continuous(labels = axis_df$CHR, breaks = axis_df$bp_cumulative) +
  labs(title = "Manhattan Plot - Rice Blast Resistance GWAS",
      subtitle = "Orange: suggestive (p<1e-4)  |  Red: genome-wide (p<1e-6)",
       x = "Chromosome", y = "-log10(p-value)") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 9))

print(manhattan_plot)
ggsave("manhattan_plot.png", manhattan_plot, width = 12, height = 5, dpi = 150)
print("Manhattan plot saved!")

# =====================================================================
# Q-Q Plot — Checks GWAS model calibration
# Expected p-values (uniform null) vs observed p-values.
# A well-calibrated GWAS follows the diagonal; inflation means
# population structure was not fully corrected.
# =====================================================================

# Compute expected -log10(p) under uniform null
n_tests <- nrow(gwas_results)
expected_log10p <- -log10(seq(1/n_tests, 1, length.out = n_tests))
observed_log10p <- sort(gwas_results$neglog10P, decreasing = TRUE)

qq_df <- data.frame(
  Expected = expected_log10p,
  Observed = observed_log10p
)

# Genomic inflation factor lambda (should be close to 1.0 = well-calibrated)
median_chisq_obs <- median(qchisq(1 - 10^(-gwas_results$neglog10P), df = 1),
                           na.rm = TRUE)
lambda <- round(median_chisq_obs / qchisq(0.5, df = 1), 3)
print(paste("Genomic inflation factor (lambda):", lambda))
# lambda close to 1.0 = good; >1.1 = inflation (population structure not fully corrected)

qq_plot <- ggplot(qq_df, aes(x = Expected, y = Observed)) +
  geom_point(size = 1.2, alpha = 0.6, color = "#2980b9") +
  geom_abline(slope = 1, intercept = 0,
              color = "red", linewidth = 0.9, linetype = "dashed") +
  geom_hline(yintercept = -log10(1e-4),
             linetype = "dotted", color = "orange", linewidth = 0.7) +
  annotate("text", x = 0.3, y = max(qq_df$Observed) * 0.95,
           label = paste0("λ = ", lambda),
           size = 5, fontface = "bold", color = "#c0392b", hjust = 0) +
  annotate("text", x = 0.3, y = max(qq_df$Observed) * 0.85,
           label = "λ ≈ 1.0 = well calibrated",
           size = 3.5, color = "gray40", hjust = 0) +
  labs(
    title    = "Q-Q Plot — GWAS of Rice Blast Resistance",
    subtitle = "Expected vs Observed -log10(p-values) | Red line = perfect calibration",
    x        = "Expected -log10(p-value)",
    y        = "Observed -log10(p-value)"
  ) +
  theme_minimal(base_size = 13)

print(qq_plot)
ggsave("qq_plot.png", qq_plot, width = 7, height = 6, dpi = 150)
print("Q-Q plot saved as qq_plot.png")

# =====================================================================
# Select Top SNPs for ML (p < 1e-4)
# =====================================================================
top_snps <- gwas_results[gwas_results$Pval < 1e-4, "SNP"]
print(paste("Top significant SNPs selected for ML:", length(top_snps)))

# Save everything
save(X_encoded, Y_ml, pca_res, pca_covs, gwas_results, gwas_plot, top_snps,
     file = "rice_blast_gwas.RData")
print("All saved - ready for ML!")