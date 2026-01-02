# ============================================================
# 04_make_figures.R
# Generate and save all figures for the project
# ============================================================

# Create figures folder if it doesn't exist
dir.create("figures", showWarnings = FALSE)

# ----------------------------
# Load required objects
# ----------------------------
beta <- readRDS("GSE42752_beta_matrix_final.rds")
y    <- readRDS("GSE42752_labels_final.rds")
cv_fit_lasso <- readRDS("cv_lasso.rds")   # if you saved it
cv_fit_en    <- readRDS("cv_enet.rds")    # if you saved it

# ----------------------------
# Figure 1: LASSO CV AUC
# ----------------------------
png("figures/lasso_cv_auc.png", width = 800, height = 600)

op <- par(no.readonly = TRUE)
par(mar = c(5, 5, 6, 2))

plot(cv_fit_lasso)
title("LASSO Logistic Regression – 5-fold CV AUC", line = 4)

par(op)
dev.off()

# ----------------------------
# Figure 2: Elastic Net CV AUC
# ----------------------------
png("figures/elastic_net_cv_auc.png", width = 800, height = 600)

op <- par(no.readonly = TRUE)
par(mar = c(5, 5, 6, 2))

plot(cv_fit_en)
title("Elastic Net Logistic Regression – 5-fold CV AUC", line = 4)

par(op)
dev.off()

# ----------------------------
# Figure 3: CpG boxplot (cg16601494)
# ----------------------------
cpg_id <- "cg16601494"

plot_df <- data.frame(
  Methylation = as.numeric(beta[cpg_id, ]),
  Tissue = factor(ifelse(y == 1, "Tumor", "Normal"),
                  levels = c("Normal", "Tumor"))
)

png("figures/cg16601494_boxplot.png", width = 800, height = 600)

boxplot(
  Methylation ~ Tissue,
  data = plot_df,
  col = c("lightblue", "salmon"),
  main = paste("DNA Methylation at", cpg_id),
  ylab = "Beta Value",
  xlab = "Tissue Type"
)

stripchart(
  Methylation ~ Tissue,
  data = plot_df,
  vertical = TRUE,
  method = "jitter",
  pch = 21,
  bg = "gray",
  add = TRUE
)

dev.off()

