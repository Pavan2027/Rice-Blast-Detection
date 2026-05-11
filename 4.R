# =====================================================================
# SCRIPT 4: ML Models — LightGBM + GBLUP + TabNet
# =====================================================================

load("rice_blast_gwas.RData")

library(lightgbm)
library(rrBLUP)
library(torch)
library(tabnet)
library(pROC)
library(ggplot2)
library(caret)

# =====================================================================
# Step 1: Prepare Feature Matrix
# =====================================================================
print(paste("Using", length(top_snps), "GWAS-selected SNPs"))

X_features <- X_encoded[, top_snps]
y_labels   <- Y_ml$Disease_Status

# Convert to numeric binary: Resistant=0, Susceptible=1
y_binary <- as.integer(y_labels == "Susceptible")

print(paste("Dataset:", nrow(X_features), "plants x", ncol(X_features), "SNPs"))
print(table(y_labels))

# =====================================================================
# Step 2: Train/Test Split (80/20)
# =====================================================================
set.seed(42)
train_idx <- createDataPartition(y_labels, p = 0.8, list = FALSE)

X_train <- X_features[train_idx, ]
X_test  <- X_features[-train_idx, ]
y_train <- y_binary[train_idx]
y_test  <- y_binary[-train_idx]
y_train_factor <- y_labels[train_idx]
y_test_factor  <- y_labels[-train_idx]

print(paste("Train:", nrow(X_train), "| Test:", nrow(X_test)))

# =====================================================================
# MODEL 1: LightGBM
# =====================================================================
print("Training LightGBM...")

# Convert to LightGBM dataset format
dtrain <- lgb.Dataset(
  data  = as.matrix(X_train),
  label = y_train
)

# Parameters
lgb_params <- list(
  objective        = "binary",
  metric           = "auc",
  learning_rate    = 0.05,
  num_leaves       = 15,
  min_data_in_leaf = 10,
  feature_fraction = 0.8,
  bagging_fraction = 0.8,
  bagging_freq     = 5,
  verbose          = -1
)

# Train with cross validation to find best number of trees
lgb_cv <- lgb.cv(
  params   = lgb_params,
  data     = dtrain,
  nfold    = 5,
  nrounds  = 500,
  early_stopping_rounds = 30,
  verbose  = -1
)

best_iter <- lgb_cv$best_iter
print(paste("LightGBM best iteration:", best_iter))

# Train final model
lgb_model <- lgb.train(
  params  = lgb_params,
  data    = dtrain,
  nrounds = best_iter,
  verbose = -1
)

# Predictions
lgb_probs <- predict(lgb_model, as.matrix(X_test))
lgb_pred  <- ifelse(lgb_probs >= 0.5, "Susceptible", "Resistant")
lgb_pred  <- factor(lgb_pred, levels = c("Resistant", "Susceptible"))

# Evaluation
lgb_cm  <- confusionMatrix(lgb_pred, y_test_factor, positive = "Susceptible")
lgb_roc <- roc(y_test_factor, lgb_probs,
               levels = c("Resistant", "Susceptible"), quiet = TRUE)

print("--- LightGBM Results ---")
print(lgb_cm)
print(paste("AUC:", round(auc(lgb_roc), 3)))

# Feature Importance
lgb_imp <- lgb.importance(lgb_model)
print("Top SNPs by LightGBM importance:")
print(head(lgb_imp, 10))

# =====================================================================
# MODEL 2: GBLUP (Genomic Best Linear Unbiased Prediction)
# =====================================================================
print("Training GBLUP...")

# GBLUP uses the full SNP matrix (not just top 19) — this is standard
# It builds a Genomic Relationship Matrix (G matrix) from all SNPs
X_gblup_full <- X_encoded  # use all 35800 SNPs for G matrix

# Build the Genomic Relationship Matrix (G matrix)
print("Building Genomic Relationship Matrix (G matrix)...")
G_matrix <- A.mat(X_gblup_full)  # from rrBLUP
print(paste("G matrix dimensions:", nrow(G_matrix), "x", ncol(G_matrix)))

# GBLUP works with continuous phenotype values
# Use raw blast scores (not binary) for prediction, then threshold
y_continuous <- Y_ml$Blast_Score

# Split G matrix same as before
G_train <- G_matrix[train_idx, train_idx]
G_test  <- G_matrix[-train_idx, train_idx]

# Fit GBLUP model using mixed.solve
gblup_fit <- mixed.solve(
  y = y_continuous[train_idx],
  K = G_train
)

# Predict test set using genomic breeding values
# u_hat = G_test_train %*% solve(G_train) %*% u_train
u_train     <- gblup_fit$u
G_train_inv <- solve(G_train + diag(1e-6, nrow(G_train)))
gblup_pred_continuous <- as.vector(G_test %*% G_train_inv %*% u_train) + 
                          as.numeric(gblup_fit$beta)

# Convert continuous predictions to binary using threshold 3.5
gblup_pred  <- factor(
  ifelse(gblup_pred_continuous >= 3.5, "Susceptible", "Resistant"),
  levels = c("Resistant", "Susceptible")
)

# Evaluation
gblup_cm  <- confusionMatrix(gblup_pred, y_test_factor, positive = "Susceptible")
gblup_roc <- roc(y_test_factor, gblup_pred_continuous,
                 levels = c("Resistant", "Susceptible"), quiet = TRUE)

print("--- GBLUP Results ---")
print(gblup_cm)
print(paste("AUC:", round(auc(gblup_roc), 3)))

# =====================================================================
# MODEL 3: TabNet
# =====================================================================
print("Training TabNet...")

# TabNet needs a dataframe with factor outcome
train_df <- as.data.frame(X_train)
test_df  <- as.data.frame(X_test)
train_df$outcome <- y_train_factor
test_df$outcome  <- y_test_factor

# Fit TabNet
set.seed(42)
tabnet_model <- tabnet_fit(
  outcome ~ .,
  data          = train_df,
  epochs        = 50,
  batch_size    = 64,
  learn_rate    = 0.02,
  num_steps     = 3,
  attention_width = 8,
  verbose       = TRUE
)

# Predictions
tabnet_probs_df <- predict(tabnet_model, test_df, type = "prob")
tabnet_probs    <- tabnet_probs_df$.pred_Susceptible
tabnet_pred     <- predict(tabnet_model, test_df, type = "class")$.pred_class

# Evaluation
tabnet_cm  <- confusionMatrix(tabnet_pred, y_test_factor, positive = "Susceptible")
tabnet_roc <- roc(y_test_factor, tabnet_probs,
                  levels = c("Resistant", "Susceptible"), quiet = TRUE)

print("--- TabNet Results ---")
print(tabnet_cm)
print(paste("AUC:", round(auc(tabnet_roc), 3)))

# =====================================================================
# Step 3: ROC Curves — All 3 Models Together
# =====================================================================
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

roc_plot <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 1.3) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray60") +
  scale_color_manual(values = c("#e74c3c", "#2980b9", "#27ae60")) +
  labs(title = "ROC Curves — Rice Blast Resistance Prediction",
       subtitle = "LightGBM vs GBLUP vs TabNet",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(roc_plot)
ggsave("roc_curves.png", roc_plot, width = 8, height = 6, dpi = 150)
print("ROC plot saved!")

# =====================================================================
# Step 4: LightGBM Feature Importance Plot
# =====================================================================
lgb_imp_plot <- ggplot(head(lgb_imp, 19),
                        aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_col(fill = "#e74c3c") +
  coord_flip() +
  labs(title = "SNP Importance — LightGBM",
       subtitle = "Ranked by information gain",
       x = "SNP", y = "Gain") +
  theme_minimal()

print(lgb_imp_plot)
ggsave("snp_importance_lgbm.png", lgb_imp_plot, width = 8, height = 6, dpi = 150)

# =====================================================================
# Step 5: Final Comparison Table
# =====================================================================
summary_table <- data.frame(
  Model = c("LightGBM", "GBLUP", "TabNet"),
  Accuracy = c(
    round(lgb_cm$overall["Accuracy"], 3),
    round(gblup_cm$overall["Accuracy"], 3),
    round(tabnet_cm$overall["Accuracy"], 3)
  ),
  AUC = c(
    round(auc(lgb_roc), 3),
    round(auc(gblup_roc), 3),
    round(auc(tabnet_roc), 3)
  ),
  Sensitivity = c(
    round(lgb_cm$byClass["Sensitivity"], 3),
    round(gblup_cm$byClass["Sensitivity"], 3),
    round(tabnet_cm$byClass["Sensitivity"], 3)
  ),
  Specificity = c(
    round(lgb_cm$byClass["Specificity"], 3),
    round(gblup_cm$byClass["Specificity"], 3),
    round(tabnet_cm$byClass["Specificity"], 3)
  )
)

print("========================================")
print("       FINAL MODEL COMPARISON           ")
print("========================================")
print(summary_table)
print("========================================")

# Save everything
save(lgb_model, gblup_fit, tabnet_model,
     lgb_cm, gblup_cm, tabnet_cm,
     lgb_roc, gblup_roc, tabnet_roc,
     summary_table,
     file = "rice_blast_final.RData")

print("Project complete! All models saved.")