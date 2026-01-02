# ============================================================
# Tumor vs Normal DNA Methylation (GSE42752)
# Minimal end-to-end pipeline:
# load -> clean CpGs -> label -> variance filter -> impute -> LASSO/EN -> plot
# ============================================================

# --- 0) Setup ---
file_path <- "C:/Users/enock_p22oyv9/OneDrive/Desktop/BIOINFORMATICIAN/BIOINFORMATICS PROJECTS/Project 3/GSE42752_series_matrix.txt.gz"

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery", update = FALSE, ask = FALSE)
if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")

library(GEOquery)
library(glmnet)

# --- 1) Load GEO series matrix from local file ---
gse  <- getGEO(filename = file_path)
eset <- if (is.list(gse)) gse[[1]] else gse

beta <- exprs(eset)        # CpGs x Samples
pheno <- pData(eset)

# --- 2) Keep only CpG probes (cg*) and ensure numeric ---
probe_ids <- gsub('"', "", rownames(beta))
keep_cpg  <- grepl("^cg", probe_ids)

beta <- beta[keep_cpg, , drop = FALSE]
rownames(beta) <- probe_ids[keep_cpg]
storage.mode(beta) <- "numeric"

# --- 3) Create Tumor/Normal labels using phenotype text (robust rules) ---
meta_text <- tolower(apply(pheno, 1, paste, collapse = " | "))

is_unrelated_normal <- grepl("cancer-unrelated", meta_text)
is_tumor  <- grepl("adenocarcinoma", meta_text)
is_normal <- (grepl("normal", meta_text) | is_unrelated_normal) & !is_tumor

y <- rep(NA_integer_, length(meta_text))
y[is_tumor]  <- 1L
y[is_normal] <- 0L
names(y) <- rownames(pheno)

keep_samples <- !is.na(y)
beta <- beta[, keep_samples, drop = FALSE]
y    <- y[keep_samples]

cat("Final beta matrix (CpGs x Samples):\n"); print(dim(beta))
cat("Class distribution (0=Normal, 1=Tumor):\n"); print(table(y))
cat("Beta range:\n"); print(range(beta, na.rm = TRUE))

# Save clean objects (optional but useful for reproducibility)
saveRDS(beta,  "GSE42752_beta_matrix_final.rds")
saveRDS(y,     "GSE42752_labels_final.rds")
saveRDS(pheno, "GSE42752_pheno.rds")

# --- 4) Variance filtering: top 10,000 CpGs ---
top_k <- 10000
cpg_var <- apply(beta, 1, var, na.rm = TRUE)
top_idx <- order(cpg_var, decreasing = TRUE)[seq_len(min(top_k, length(cpg_var)))]

beta_top <- beta[top_idx, , drop = FALSE]
cat("Variance-filtered beta (CpGs x Samples):\n"); print(dim(beta_top))
saveRDS(beta_top, "GSE42752_beta_matrix_top10k.rds")

# --- 5) Build modeling matrix: samples x CpGs + fix names + impute ---
X <- t(beta_top)                 # Samples x CpGs
X <- as.matrix(X)
storage.mode(X) <- "numeric"

# Fix any bad/duplicated feature names (rare but safe)
bad_names <- is.na(colnames(X)) | colnames(X) == ""
if (any(bad_names)) colnames(X)[bad_names] <- paste0("CpG_", which(bad_names))
colnames(X) <- make.unique(colnames(X))

# Replace any non-finite with NA, then CpG-wise mean imputation
X[is.nan(X)] <- NA
X[is.infinite(X)] <- NA

na_cols <- which(colSums(is.na(X)) > 0)
if (length(na_cols) > 0) {
  for (j in na_cols) {
    X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)
  }
}

cat("Validation (must be 0): NA=", sum(is.na(X)),
    " NaN=", sum(is.nan(X)),
    " Inf=", sum(is.infinite(X)), "\n")

# --- 6) LASSO logistic regression (5-fold CV AUC) ---
set.seed(123)
cv_lasso <- cv.glmnet(
  x = X, y = y,
  family = "binomial",
  alpha = 1,
  nfolds = 5,
  type.measure = "auc",
  standardize = TRUE
)

cat("\nLASSO CV AUC (best):\n"); print(max(cv_lasso$cvm))
cat("LASSO lambda.min:\n"); print(cv_lasso$lambda.min)
cat("LASSO lambda.1se:\n"); print(cv_lasso$lambda.1se)

coef_lasso <- coef(cv_lasso, s = "lambda.1se")
sel_lasso <- rownames(coef_lasso)[coef_lasso[,1] != 0]
sel_lasso <- setdiff(sel_lasso, "(Intercept)")

cat("Selected CpGs (LASSO, lambda.1se):", length(sel_lasso), "\n")
print(head(sel_lasso, 20))

# Fix title overlap for LASSO CV plot
op <- par(no.readonly = TRUE)

par(mar = c(5, 5, 7, 2))   # increase TOP margin more

plot(cv_fit)

title(
  "LASSO Logistic Regression – 5-fold CV AUC",
  line = 5,               # push title higher
  cex.main = 1.2
)

par(op)


# --- 7) Elastic Net logistic regression (alpha = 0.5) ---
set.seed(123)
cv_enet <- cv.glmnet(
  x = X, y = y,
  family = "binomial",
  alpha = 0.5,
  nfolds = 5,
  type.measure = "auc",
  standardize = TRUE
)

cat("\nElastic Net CV AUC (best):\n"); print(max(cv_enet$cvm))
cat("Elastic Net lambda.min:\n"); print(cv_enet$lambda.min)
cat("Elastic Net lambda.1se:\n"); print(cv_enet$lambda.1se)

coef_enet <- coef(cv_enet, s = "lambda.1se")
sel_enet <- rownames(coef_enet)[coef_enet[,1] != 0]
sel_enet <- setdiff(sel_enet, "(Intercept)")

cat("Selected CpGs (Elastic Net, lambda.1se):", length(sel_enet), "\n")
print(head(sel_enet, 20))


op <- par(no.readonly = TRUE)

par(mar = c(5, 5, 6, 2))  # increase top margin

plot(cv_fit_en)

title(
  "Elastic Net Logistic Regression – 5-fold CV AUC",
  line = 4,
  cex.main = 1.2
)

par(op)

# --- 8) Plot methylation for a selected CpG (use cg16601494 if present) ---
cpg_id <- if ("cg16601494" %in% rownames(beta)) "cg16601494" else sel_lasso[1]

plot_df <- data.frame(
  Methylation = as.numeric(beta[cpg_id, ]),
  Tissue = factor(ifelse(y == 1, "Tumor", "Normal"), levels = c("Normal", "Tumor"))
)

boxplot(
  Methylation ~ Tissue, data = plot_df,
  col = c("lightblue", "salmon"),
  main = paste("DNA Methylation at", cpg_id),
  ylab = "Beta Value", xlab = "Tissue Type"
)
stripchart(
  Methylation ~ Tissue, data = plot_df,
  vertical = TRUE, method = "jitter",
  pch = 21, bg = "gray", add = TRUE
)


