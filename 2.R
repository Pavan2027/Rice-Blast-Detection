# =====================================================================
# SCRIPT 2: SNP Encoding — Convert Letters to Numbers
# =====================================================================


# Load Script 1 output — instant, no waiting
load("rice_blast_ready.RData")

# Confirm it loaded
print(paste("Plants:", nrow(X_ml)))
print(paste("SNPs:", ncol(X_ml)))

print("Encoding SNPs... this may take 1-2 minutes for 36,901 markers")

# Each column is one SNP — find the most common allele (reference)
# then encode: ref/ref = 0, mixed = 1, alt/alt = 2

encode_snp_column <- function(snp_col) {
  # Get unique non-NA values
  alleles <- na.omit(unique(snp_col))
  
  # Find reference allele (most frequent)
  freq_table <- table(snp_col)
  ref_allele <- names(which.max(freq_table))
  
  # Encode
  encoded <- ifelse(snp_col == ref_allele, 0L,
             ifelse(is.na(snp_col), NA_integer_, 2L))
  return(encoded)
}

# Apply to all 36,901 SNP columns
X_encoded <- apply(X_ml, 2, encode_snp_column)
rownames(X_encoded) <- rownames(X_ml)

print(paste("Encoding done. Dimensions:", dim(X_encoded)[1], "x", dim(X_encoded)[2]))

# Verify — should now show 0s and 2s
print("Sample of encoded values:")
print(X_encoded[1:5, 1:5])

# =====================================================================
# Quality Control: Remove low-quality SNPs
# =====================================================================

# Filter 1: Remove SNPs with >10% missing data
miss_rate <- apply(X_encoded, 2, function(x) mean(is.na(x)))
X_encoded <- X_encoded[, miss_rate <= 0.10]
print(paste("SNPs after missingness filter:", ncol(X_encoded)))

# Filter 2: Remove SNPs with Minor Allele Frequency < 5%
maf <- apply(X_encoded, 2, function(x) {
  freq <- mean(x, na.rm = TRUE) / 2
  min(freq, 1 - freq)
})
X_encoded <- X_encoded[, maf >= 0.05]
print(paste("SNPs after MAF filter:", ncol(X_encoded)))

# Filter 3: Impute remaining NAs with column mean
for (j in seq_len(ncol(X_encoded))) {
  na_idx <- is.na(X_encoded[, j])
  if (any(na_idx)) {
    X_encoded[na_idx, j] <- round(mean(X_encoded[, j], na.rm = TRUE))
  }
}

print("No more NAs after imputation.")
print(paste("Final dataset:", nrow(X_encoded), "plants x", ncol(X_encoded), "SNPs"))

# =====================================================================
# Save the encoded clean version
# =====================================================================
save(X_encoded, Y_ml, file = "rice_blast_encoded.RData")
print("Saved as rice_blast_encoded.RData — use this for all ML steps")