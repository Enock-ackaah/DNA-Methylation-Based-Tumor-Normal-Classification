# DNA-Methylation-Based-Tumor-Normal-Classification
Interpretable DNA methylation–based models for tumor vs normal classification in colon cancer using LASSO and Elastic Net logistic regression.
Methods
Study Design

This study evaluates whether genome-wide DNA methylation profiles can distinguish colon adenocarcinoma tissue from normal colon tissue. A supervised classification framework was adopted, emphasizing interpretability, robustness, and suitability for high-dimensional epigenomic data.

Data Source

DNA methylation data were obtained from the Gene Expression Omnibus (GEO) under accession GSE42752. The dataset consists of colon adenocarcinoma tissue samples and corresponding normal colon tissue samples profiled using the Illumina HumanMethylation450 BeadChip, which measures methylation beta values at approximately 450,000 CpG sites genome-wide.

Data Loading and Preprocessing

The GEO series matrix file was loaded locally using the GEOquery package. Expression values were extracted as beta values (ranging from 0 to 1), where higher values indicate greater methylation.

Preprocessing steps included:

Removal of non-CpG probes by retaining probes with identifiers beginning with “cg”

Enforcement of numeric data types

Validation of matrix dimensions and value ranges

Extraction and parsing of phenotype metadata to identify tumor and normal tissue samples

Tumor versus normal labels were assigned using keyword-based parsing of GEO sample annotations, with careful exclusion of ambiguous or mislabeled samples.

Feature Reduction

Given the high-dimensional nature of DNA methylation data (p ≫ n), variance-based feature filtering was applied prior to modeling. CpG-wise variance was computed across samples, and the top 10,000 most variable CpG sites were retained. This approach reduces noise, improves numerical stability, and preserves features most likely to reflect biologically meaningful differences between tumor and normal tissue.

Handling Missing and Invalid Values

DNA methylation data may contain missing or invalid values due to technical limitations. The following steps were applied:

Removal of CpG features with zero variance

Replacement of non-finite values (NA, NaN, Inf) using CpG-wise mean imputation

Verification that the final design matrix contained no missing or invalid values prior to modeling

These steps ensured compatibility with regularized regression algorithms.

Statistical Modeling
LASSO Logistic Regression

A LASSO (Least Absolute Shrinkage and Selection Operator) logistic regression model was fitted to classify tumor versus normal tissue. LASSO applies an L1 penalty to regression coefficients, enforcing sparsity and enabling automatic feature selection.

Key modeling details:

Binary outcome: tumor = 1, normal = 0

Five-fold cross-validation

Performance metric: Area Under the ROC Curve (AUC)

Conservative model selection using the λ₁SE criterion

Elastic Net Logistic Regression

To assess robustness and account for correlated CpG features, an Elastic Net logistic regression model (α = 0.5) was also fitted. Elastic Net combines LASSO and Ridge penalties, allowing groups of correlated features to be selected together while still controlling overfitting.

Model Evaluation

Model performance was evaluated using five-fold cross-validation. The primary evaluation metric was the area under the receiver operating characteristic curve (AUC). Stability was assessed by examining performance across a wide range of regularization strengths and by comparing feature selection results between LASSO and Elastic Net models.

Visualization

For interpretability, methylation beta values at selected CpG sites were visualized using boxplots stratified by tissue type. These plots provided intuitive confirmation of the statistical findings and highlighted the biological separation between tumor and normal samples.

Software and Tools

All analyses were conducted in R using the following packages:

GEOquery (data retrieval and parsing)

glmnet (regularized logistic regression)

Base R functions for preprocessing, visualization, and validation

Reproducibility

All preprocessing, modeling, and evaluation steps are fully scripted and reproducible. Intermediate datasets (cleaned matrices and labels) are saved as RDS files to ensure consistency across analysis stages.

Ethical Considerations
