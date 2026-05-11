# =====================================================================
# SCRIPT 5: Final Outputs — EDA + ML + Metrics + Resistance Scale
# =====================================================================


load("rice_blast_gwas.RData")
load("rice_blast_final.RData")

library(ggplot2)
library(dplyr)
library(pROC)
library(caret)
library(lightgbm)

# =====================================================================
# PART 1: EDA — Exploratory Data Analysis
# =====================================================================

# ------- Plot 1: Blast Score Distribution -------
p1 <- ggplot(Y_ml, aes(x = Blast_Score)) +
  geom_histogram(aes(fill = Disease_Status), bins = 10,
                 color = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("Resistant" = "#27ae60",
                                "Susceptible" = "#e74c3c")) +
  geom_vline(xintercept = 3.5, linetype = "dashed",
             color = "black", linewidth = 0.8) +
  annotate("text", x = 1.8, y = 45, label = "Resistant\n(Score ≤ 3)",
           color = "#27ae60", fontface = "bold", size = 4) +
  annotate("text", x = 5.5, y = 45, label = "Susceptible\n(Score ≥ 4)",
           color = "#e74c3c", fontface = "bold", size = 4) +
  labs(title = "Distribution of Rice Blast Disease Scores",
       subtitle = "Shanghang Field Trial | 364 Rice Accessions",
       x = "Blast Score (0-9)", y = "Number of Plants",
       fill = "Disease Status") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

print(p1)
ggsave("eda_1_blast_distribution.png", p1, width = 8, height = 5, dpi = 150)

# ------- Plot 2: Class Balance -------
class_df <- as.data.frame(table(Y_ml$Disease_Status))
colnames(class_df) <- c("Status", "Count")
class_df$Percentage <- round(class_df$Count / sum(class_df$Count) * 100, 1)
class_df$Label <- paste0(class_df$Count, "\n(", class_df$Percentage, "%)")

p2 <- ggplot(class_df, aes(x = Status, y = Count, fill = Status)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_text(aes(label = Label), vjust = -0.5, fontface = "bold", size = 5) +
  scale_fill_manual(values = c("Resistant" = "#27ae60",
                                "Susceptible" = "#e74c3c")) +
  scale_y_continuous(limits = c(0, 240)) +
  labs(title = "Class Balance — Resistant vs Susceptible",
       subtitle = "Nearly balanced dataset — ideal for ML",
       x = "", y = "Number of Accessions") +
  theme_minimal(base_size = 13)

print(p2)
ggsave("eda_2_class_balance.png", p2, width = 6, height = 5, dpi = 150)

# ------- Plot 3: Blast Score by Subpopulation -------
subpop_df <- Y_ml %>%
  select(Disease_Status, Blast_Score, Subpopulations) %>%
  filter(Subpopulations != "Not available")

p3 <- ggplot(subpop_df, aes(x = reorder(Subpopulations, Blast_Score, median),
                              y = Blast_Score, fill = Subpopulations)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = 21,
               outlier.size = 1.5, alpha = 0.85) +
  geom_hline(yintercept = 3.5, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  annotate("text", x = 0.6, y = 4.2, label = "Threshold",
           color = "red", size = 3.5) +
  labs(title = "Blast Resistance Score by Rice Subpopulation",
       subtitle = "Lower score = more resistant",
       x = "Subpopulation", y = "Blast Score") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p3)
ggsave("eda_3_subpopulation_blast.png", p3, width = 8, height = 5, dpi = 150)

# ------- Plot 4: PCA colored by Subpopulation -------
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$SampleID <- rownames(pca_df)
pca_df$Disease_Status <- Y_ml$Disease_Status
pca_df$Subpopulation  <- Y_ml$Subpopulations

variance_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

p4 <- ggplot(pca_df %>% filter(Subpopulation != "Not available"),
             aes(x = PC1, y = PC2, color = Subpopulation,
                 shape = Disease_Status)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_shape_manual(values = c("Resistant" = 16, "Susceptible" = 4)) +
  labs(title = "PCA — Population Structure of Rice Panel",
       subtitle = "Shape = Disease Status | Color = Subpopulation",
       x = paste0("PC1 (", variance_exp[1], "% variance)"),
       y = paste0("PC2 (", variance_exp[2], "% variance)"),
       color = "Subpopulation", shape = "Disease Status") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

print(p4)
ggsave("eda_4_pca_subpopulation.png", p4, width = 9, height = 6, dpi = 150)

# ------- Plot 5: Top GWAS SNPs Bar Chart -------
top10_gwas <- head(gwas_results[order(gwas_results$Pval), ], 10)
top10_gwas$SNP <- factor(top10_gwas$SNP,
                          levels = top10_gwas$SNP[order(top10_gwas$neglog10P)])

p5 <- ggplot(top10_gwas, aes(x = SNP, y = neglog10P, fill = neglog10P)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = -log10(1e-4), linetype = "dashed",
             color = "orange", linewidth = 0.7) +
  scale_fill_gradient(low = "#f39c12", high = "#c0392b") +
  coord_flip() +
  labs(title = "Top 10 GWAS SNPs — Rice Blast Resistance",
       subtitle = "Orange line = suggestive significance threshold (p < 1e-4)",
       x = "SNP ID", y = "-log10(p-value)") +
  theme_minimal(base_size = 12)

print(p5)
ggsave("eda_5_top_gwas_snps.png", p5, width = 8, height = 5, dpi = 150)

print("EDA plots done — 5 plots saved!")

# =====================================================================
# PART 2: ML MODEL — Re-prepare data and run predictions
# =====================================================================

# Rebuild train/test split (same seed = same split as Script 4)
X_features <- X_encoded[, top_snps]
y_labels   <- Y_ml$Disease_Status
y_binary   <- as.integer(y_labels == "Susceptible")

set.seed(42)
train_idx      <- createDataPartition(y_labels, p = 0.8, list = FALSE)
X_test         <- X_features[-train_idx, ]
y_test_factor  <- y_labels[-train_idx]
y_test_binary  <- y_binary[-train_idx]

# Get predictions from saved models
lgb_probs    <- predict(lgb_model, as.matrix(X_test))
lgb_pred     <- factor(ifelse(lgb_probs >= 0.5, "Susceptible", "Resistant"),
                        levels = c("Resistant", "Susceptible"))

# GBLUP predictions (rebuild from saved fit)
G_matrix  <- rrBLUP::A.mat(X_encoded)
G_test    <- G_matrix[-train_idx, train_idx]
G_train   <- G_matrix[train_idx, train_idx]
G_train_inv <- solve(G_train + diag(1e-6, nrow(G_train)))
gblup_pred_continuous <- as.vector(G_test %*% G_train_inv %*% gblup_fit$u) +
                          as.numeric(gblup_fit$beta)
gblup_pred <- factor(ifelse(gblup_pred_continuous >= 3.5,
                              "Susceptible", "Resistant"),
                      levels = c("Resistant", "Susceptible"))

# TabNet predictions
test_df          <- as.data.frame(X_test)
test_df$outcome  <- y_test_factor
tabnet_probs_df  <- predict(tabnet_model, test_df, type = "prob")
tabnet_probs     <- tabnet_probs_df$.pred_Susceptible
tabnet_pred      <- predict(tabnet_model, test_df, type = "class")$.pred_class

# =====================================================================
# PART 3: PERFORMANCE METRICS
# =====================================================================

# Confusion matrices
lgb_cm    <- confusionMatrix(lgb_pred,    y_test_factor, positive = "Susceptible")
gblup_cm  <- confusionMatrix(gblup_pred,  y_test_factor, positive = "Susceptible")
tabnet_cm <- confusionMatrix(tabnet_pred, y_test_factor, positive = "Susceptible")

# ROC + AUC
lgb_roc    <- roc(y_test_factor, lgb_probs,
                   levels = c("Resistant","Susceptible"), quiet = TRUE)
gblup_roc  <- roc(y_test_factor, gblup_pred_continuous,
                   levels = c("Resistant","Susceptible"), quiet = TRUE)
tabnet_roc <- roc(y_test_factor, tabnet_probs,
                   levels = c("Resistant","Susceptible"), quiet = TRUE)

# ------- Final Metrics Table -------
metrics_table <- data.frame(
  Model = c("LightGBM", "GBLUP", "TabNet"),
  Accuracy    = c(round(lgb_cm$overall["Accuracy"], 3),
                  round(gblup_cm$overall["Accuracy"], 3),
                  round(tabnet_cm$overall["Accuracy"], 3)),
  AUC         = c(round(auc(lgb_roc), 3),
                  round(auc(gblup_roc), 3),
                  round(auc(tabnet_roc), 3)),
  Sensitivity = c(round(lgb_cm$byClass["Sensitivity"], 3),
                  round(gblup_cm$byClass["Sensitivity"], 3),
                  round(tabnet_cm$byClass["Sensitivity"], 3)),
  Specificity = c(round(lgb_cm$byClass["Specificity"], 3),
                  round(gblup_cm$byClass["Specificity"], 3),
                  round(tabnet_cm$byClass["Specificity"], 3)),
  F1_Score    = c(round(lgb_cm$byClass["F1"], 3),
                  round(gblup_cm$byClass["F1"], 3),
                  round(tabnet_cm$byClass["F1"], 3)),
  Kappa       = c(round(lgb_cm$overall["Kappa"], 3),
                  round(gblup_cm$overall["Kappa"], 3),
                  round(tabnet_cm$overall["Kappa"], 3))
)

print("==========================================")
print("       FINAL MODEL COMPARISON TABLE       ")
print("==========================================")
print(metrics_table)
print("==========================================")

# ------- Plot 6: Metrics Comparison Bar Chart -------
library(tidyr)
metrics_long <- metrics_table %>%
  pivot_longer(cols = -Model,
               names_to = "Metric", values_to = "Value")

p6 <- ggplot(metrics_long, aes(x = Metric, y = Value, fill = Model)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_text(aes(label = Value), position = position_dodge(width = 0.65),
            vjust = -0.4, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("LightGBM" = "#e74c3c",
                                "GBLUP"    = "#2980b9",
                                "TabNet"   = "#27ae60")) +
  scale_y_continuous(limits = c(0, 1.1)) +
  labs(title = "Model Performance Comparison",
       subtitle = "LightGBM vs GBLUP vs TabNet — Rice Blast Resistance",
       x = "Metric", y = "Score", fill = "Model") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 20, hjust = 1))

print(p6)
ggsave("metrics_comparison.png", p6, width = 9, height = 6, dpi = 150)

# ------- Plot 7: ROC Curves -------
roc_df <- rbind(
  data.frame(FPR   = 1 - lgb_roc$specificities,
             TPR   = lgb_roc$sensitivities,
             Model = paste0("LightGBM (AUC=", round(auc(lgb_roc), 3), ")")),
  data.frame(FPR   = 1 - gblup_roc$specificities,
             TPR   = gblup_roc$sensitivities,
             Model = paste0("GBLUP (AUC=", round(auc(gblup_roc), 3), ")")),
  data.frame(FPR   = 1 - tabnet_roc$specificities,
             TPR   = tabnet_roc$sensitivities,
             Model = paste0("TabNet (AUC=", round(auc(tabnet_roc), 3), ")"))
)

p7 <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 1.3) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray60") +
  scale_color_manual(values = c("#e74c3c", "#2980b9", "#27ae60")) +
  labs(title = "ROC Curves — Rice Blast Resistance Prediction",
       subtitle = "LightGBM vs GBLUP vs TabNet",
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(p7)
ggsave("roc_curves_final.png", p7, width = 8, height = 6, dpi = 150)

# =====================================================================
# PART 4: RESISTANCE SCALE
# =====================================================================
# Use average predicted probability across all 3 models
# as a unified "Blast Resistance Score" from 0 to 100

# Get probabilities for ALL 364 plants (train + test combined)
all_lgb_probs <- predict(lgb_model, as.matrix(X_features))

all_test_df          <- as.data.frame(X_features)
all_test_df$outcome  <- y_labels
all_tabnet_probs     <- predict(tabnet_model, all_test_df,
                                type = "prob")$.pred_Susceptible

G_all     <- G_matrix
G_all_inv <- solve(G_matrix[train_idx, train_idx] +
                    diag(1e-6, length(train_idx)))
all_gblup_probs_raw <- as.vector(
  G_matrix[, train_idx] %*% G_all_inv %*% gblup_fit$u
) + as.numeric(gblup_fit$beta)

# Normalize GBLUP to 0-1 scale
all_gblup_probs <- (all_gblup_probs_raw - min(all_gblup_probs_raw)) /
                    (max(all_gblup_probs_raw) - min(all_gblup_probs_raw))

# Ensemble: average of all 3 models
ensemble_prob <- (all_lgb_probs + all_gblup_probs + all_tabnet_probs) / 3

# Resistance Score = inverse of susceptibility probability (0-100)
resistance_score <- round((1 - ensemble_prob) * 100, 1)

# Build final resistance table
resistance_df <- data.frame(
  Plant_ID       = rownames(X_features),
  Blast_Score    = Y_ml$Blast_Score,
  True_Status    = Y_ml$Disease_Status,
  Susceptibility_Prob = round(ensemble_prob, 3),
  Resistance_Score    = resistance_score,
  Resistance_Class    = cut(resistance_score,
    breaks = c(0, 20, 40, 60, 80, 100),
    labels = c("Highly Susceptible",
               "Moderately Susceptible",
               "Intermediate",
               "Moderately Resistant",
               "Highly Resistant"),
    include.lowest = TRUE)
)

# Sort by most resistant first
resistance_df <- resistance_df[order(resistance_df$Resistance_Score,
                                      decreasing = TRUE), ]

print("==========================================")
print("       TOP 15 MOST RESISTANT PLANTS       ")
print("==========================================")
print(head(resistance_df, 15))

print("==========================================")
print("       TOP 15 MOST SUSCEPTIBLE PLANTS     ")
print("==========================================")
print(tail(resistance_df, 15))

print("Resistance Class Distribution:")
print(table(resistance_df$Resistance_Class))

# ------- Plot 8: Resistance Scale Distribution -------
p8 <- ggplot(resistance_df, aes(x = Resistance_Score,
                                 fill = Resistance_Class)) +
  geom_histogram(bins = 25, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "Highly Susceptible"     = "#c0392b",
    "Moderately Susceptible" = "#e67e22",
    "Intermediate"           = "#f1c40f",
    "Moderately Resistant"   = "#2ecc71",
    "Highly Resistant"       = "#1a7a4a"
  )) +
  labs(title = "Blast Resistance Score Distribution",
       subtitle = "Ensemble Score from LightGBM + GBLUP + TabNet | 364 Rice Accessions",
       x = "Resistance Score (0 = Susceptible → 100 = Resistant)",
       y = "Number of Plants",
       fill = "Resistance Class") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

print(p8)
ggsave("resistance_scale.png", p8, width = 10, height = 6, dpi = 150)

# ------- Plot 9: Resistance Score vs Raw Blast Score -------
p9 <- ggplot(resistance_df,
             aes(x = Blast_Score, y = Resistance_Score,
                 color = Resistance_Class)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", color = "gray40",
              linetype = "dashed", se = FALSE) +
  scale_color_manual(values = c(
    "Highly Susceptible"     = "#c0392b",
    "Moderately Susceptible" = "#e67e22",
    "Intermediate"           = "#f1c40f",
    "Moderately Resistant"   = "#2ecc71",
    "Highly Resistant"       = "#1a7a4a"
  )) +
  labs(title = "Model Resistance Score vs Raw Blast Score",
       subtitle = "Validates that our ensemble score matches field observations",
       x = "Raw Field Blast Score (0-9)",
       y = "Ensemble Resistance Score (0-100)",
       color = "Resistance Class") +
  theme_minimal(base_size = 12)

print(p9)
ggsave("resistance_vs_raw_score.png", p9, width = 9, height = 6, dpi = 150)

# ------- Save Full Resistance Table as CSV -------
write.csv(resistance_df,
          "rice_blast_resistance_rankings.csv",
          row.names = FALSE)

print("==========================================")
print("ALL OUTPUTS SAVED:")
print("  eda_1_blast_distribution.png")
print("  eda_2_class_balance.png")
print("  eda_3_subpopulation_blast.png")
print("  eda_4_pca_subpopulation.png")
print("  eda_5_top_gwas_snps.png")
print("  metrics_comparison.png")
print("  roc_curves_final.png")
print("  resistance_scale.png")
print("  resistance_vs_raw_score.png")
print("  rice_blast_resistance_rankings.csv")
print("==========================================")
print("PROJECT COMPLETE!")