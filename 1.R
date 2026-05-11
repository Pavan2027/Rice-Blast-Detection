# =====================================================================
# SCRIPT 1: Data Merging for Rice Blast Disease Prediction
# =====================================================================

# Make file paths robust when running from any folder.
script_path <- NULL
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
} else {
  ofile <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) {
    script_path <- ofile
  }
}

if (!is.null(script_path) && nzchar(script_path)) {
  script_path <- gsub("~\\+~", " ", script_path)
  setwd(dirname(normalizePath(script_path)))
}

print(paste("Working directory:", getwd()))
list.files()

# =====================================================================
# Step 1: Load Both Datasets
# =====================================================================
#install.packages("readxl")  # Skip if already installed
library(readxl)

print("Loading datasets... genotype file will take 15-30 seconds")

# 1A. Load Disease Scores (Y) — FIXED: read_excel not read.excel
# Load from Sheet 2
blast_data <- read_excel("data/12284_2016_116_MOESM1_ESM.xlsx", sheet = 2, skip = 2)
blast_data <- as.data.frame(blast_data)

# Verify it loaded correctly
class(blast_data)   # Must say "data.frame"
dim(blast_data)     # Must show rows × columns
head(blast_data)    # Must show actual data
colnames(blast_data)  # Read this carefully — find the NSFTV ID column name

# 1B. Load 44K DNA Markers (X)
genotypes <- read.csv("data/RiceDiversity.44K.MSU6.Genotypes.csv.gz", 
                      header = TRUE, row.names = 1)

# Quick size check
dim(blast_data)
dim(genotypes)

# =====================================================================
# Step 2: Set Blast Column Name (Exact name from your file)
# =====================================================================
BLAST_COLUMN_NAME <- "Shanghang (Fujian province, Southeast of China)"

# =====================================================================
# Step 3: Clean Phenotype Data
# =====================================================================

# Remove placeholder rows
blast_data <- blast_data[blast_data$NSFTV.ID != "Not available", ]

# Trim hidden spaces
blast_data$NSFTV.ID <- trimws(as.character(blast_data$NSFTV.ID))
blast_data$NSFTV.ID <- paste0("NSFTV_", blast_data$NSFTV.ID)

# Extract blast score
blast_data$Blast_Score <- as.numeric(blast_data[[BLAST_COLUMN_NAME]])

# Drop unscored plants
blast_data <- blast_data[!is.na(blast_data$Blast_Score), ]

# Create ML Target: 0-3 = Resistant, 4-9 = Susceptible
blast_data$Disease_Status <- as.factor(
  ifelse(blast_data$Blast_Score <= 3, "Resistant", "Susceptible")
)

print("Phenotype data cleaned.")
print(table(blast_data$Disease_Status))

# =====================================================================
# Step 4: Prepare Genotype Data
# =====================================================================
geno_flipped <- t(genotypes)
rownames(geno_flipped) <- trimws(as.character(rownames(geno_flipped)))
geno_flipped <- geno_flipped[!(rownames(geno_flipped) %in% c("chr", "position")), , drop = FALSE]
print("Genotype data transposed.")

# =====================================================================
# Step 5: Merge
# =====================================================================
common_plants <- intersect(blast_data$NSFTV.ID, rownames(geno_flipped))
print(paste("Matched plants:", length(common_plants)))

if (length(common_plants) == 0) {
  stop("No matched plants found after ID normalization. Check NSFTV.ID format and genotype rownames.")
}

Y_ml <- blast_data[match(common_plants, blast_data$NSFTV.ID), ]
X_ml <- geno_flipped[common_plants, ]

# =====================================================================
# Step 6: Final Verification
# =====================================================================
print("--- Data Pipeline Complete ---")
print(paste("X_ml:", dim(X_ml)[1], "Plants,", dim(X_ml)[2], "DNA Markers"))
print(paste("Y_ml:", dim(Y_ml)[1], "Plants,", dim(Y_ml)[2], "Columns"))
print(paste("CRITICAL ALIGNMENT CHECK ->", all(rownames(X_ml) == Y_ml$NSFTV.ID)))

print("Class Distribution:")
print(table(Y_ml$Disease_Status))
print(round(prop.table(table(Y_ml$Disease_Status)) * 100, 1))

print("Sample genotype values:")
print(X_ml[1:5, 1:5])

# =====================================================================
# Step 7: Save
# =====================================================================
save(X_ml, Y_ml, file = "rice_blast_ready.RData")
print("Saved! Next time just run: load('rice_blast_ready.RData')")
