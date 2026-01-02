# Tumor vs Normal Classification Using DNA Methylation

This project investigates whether genome-wide DNA methylation profiles can distinguish colon adenocarcinoma tissue from normal colon tissue using interpretable statistical and machine-learning models.

Using publicly available Illumina HumanMethylation450 data, DNA methylation beta values were analyzed and regularized logistic regression models were trained to identify robust CpG-level epigenetic biomarkers associated with tumor status.

---

## Dataset

- **Source:** Gene Expression Omnibus (GEO)
- **Accession:** GSE42752
- **Organism:** *Homo sapiens*
- **Platform:** Illumina HumanMethylation450 BeadChip
- **Samples:** 63 colon tissue samples
- **Labels:**
  - `1` → Tumor (colon adenocarcinoma)
  - `0` → Normal (including cancer-unrelated normal tissue)
- **Input:** DNA methylation beta values (range: 0–1) at ~450,000 CpG sites

---

## Study Design

This study adopts a supervised classification framework to evaluate whether DNA methylation patterns can reliably separate tumor from normal colon tissue. The analysis emphasizes interpretability, robustness, and suitability for high-dimensional epigenomic data.

---

## Data Loading and Preprocessing

- Loaded GEO series matrix files using the **GEOquery** package
- Retained only CpG probes with identifiers beginning with `"cg"`
- Enforced numeric data types and validated beta value ranges
- Parsed phenotype metadata to assign tumor vs normal labels
- Excluded ambiguous or mislabeled samples

---

## Feature Reduction

Given the high-dimensional nature of DNA methylation data (p ≫ n), variance-based feature filtering was applied prior to modeling.

- CpG-wise variance was computed across samples
- The **top 10,000 most variable CpG sites** were retained
- This step reduces noise and improves numerical stability while preserving biologically informative features

---

## Handling Missing and Invalid Values

- Removed CpG features with zero variance
- Replaced non-finite values (NA, NaN, Inf) using CpG-wise mean imputation
- Verified that the final design matrix contained no missing or invalid values prior to modeling

---

## Methods

### LASSO Logistic Regression

- Regularized logistic regression with L1 penalty
- Automatic feature selection
- Five-fold cross-validation
- Performance metric: Area Under the ROC Curve (AUC)
- Conservative model selection using the **λ₁SE** criterion

### Elastic Net Logistic Regression

- Combined L1 and L2 penalties (α = 0.5)
- Allows correlated CpG features to be selected together
- Used to assess robustness of feature selection and model performance

---

## Results

- **LASSO Logistic Regression**
  - Cross-validated AUC: **1.0**
  - Selected **1 CpG site** (cg16601494)

- **Elastic Net Logistic Regression**
  - Cross-validated AUC: **1.0**
  - Selected **8 CpG sites**
  - cg16601494 consistently retained

Cross-validation curves demonstrated stable near-perfect performance across a wide range of regularization strengths, indicating that classification accuracy is not dependent on model complexity.

---

## Interpretation

- A single CpG site (cg16601494) shows a strong and non-overlapping methylation difference between tumor and normal tissue
- Consistent feature selection across LASSO and Elastic Net demonstrates robustness
- Results reflect a genuine biological signal rather than overfitting

---

## Visualization

- Boxplots of methylation beta values for selected CpG sites
- Clear separation between tumor and normal tissue
- Visual confirmation of statistical findings

---

## Software and Tools

- **R**
- **Bioconductor**
  - GEOquery
- **glmnet**
- Base R functions for preprocessing, modeling, and visualization

---

## Reproducibility

All preprocessing, modeling, and evaluation steps are fully scripted and reproducible. Intermediate datasets (cleaned methylation matrices and outcome labels) are saved as RDS files to ensure consistency across analysis stages.

---

## Ethical Considerations

All data used in this project are publicly available and fully de-identified. No human subjects were directly involved in this research.
