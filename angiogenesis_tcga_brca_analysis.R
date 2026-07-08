# Load packages

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(edgeR)
library(msigdbr)
library(survival)
library(survminer)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(splines)

# Download TCGA-BRCA

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)

brca_data <- GDCprepare(query)

# Extract expression matrix

counts <- assay(brca_data)

gene_symbols <- rowData(brca_data)$gene_name

rownames(counts) <- gene_symbols

counts <- counts[
  !is.na(rownames(counts)) &
    rownames(counts) != "",
]

counts <- counts[
  !duplicated(rownames(counts)),
]

# Normalize expression

dge <- DGEList(counts)

keep <- rowSums(cpm(dge) > 1) >= 10

dge <- dge[keep, ]

dge <- calcNormFactors(dge)

logCPM <- cpm(
  dge,
  log = TRUE,
  prior.count = 1
)

# Hallmark angiogenesis genes

angiogenesis_genes <- msigdbr(
  species = "Homo sapiens",
  category = "H"
) %>%
  filter(
    gs_name == "HALLMARK_ANGIOGENESIS"
  ) %>%
  pull(gene_symbol) %>%
  unique()

common_genes <- intersect(
  angiogenesis_genes,
  rownames(logCPM)
)

angiogenesis_matrix <- logCPM[
  common_genes,
]

# Clinical data

clinical <- as.data.frame(
  colData(brca_data)
)

clinical$time <- ifelse(
  is.na(clinical$days_to_death),
  clinical$days_to_last_follow_up,
  clinical$days_to_death
)

clinical$status <- ifelse(
  clinical$vital_status == "Dead",
  1,
  0
)

clinical$time_years <- clinical$time / 365.25

clinical$subtype <- factor(
  clinical$paper_BRCA_Subtype_PAM50
)

clinical$age_z <- as.numeric(
  scale(clinical$age_at_diagnosis)
)

# Align expression and clinical data

common_ids <- intersect(
  colnames(angiogenesis_matrix),
  rownames(clinical)
)

expr <- angiogenesis_matrix[
  ,
  common_ids
]

clinical <- clinical[
  common_ids,
]

# Correlation heatmap

expr_z <- t(
  scale(
    t(expr)
  )
)

pheatmap(
  cor(expr_z),
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Gene-Gene Correlation"
)

# Expression heatmap

pheatmap(
  expr_z,
  show_colnames = FALSE,
  main = "Angiogenesis Gene Expression"
)

# PCA

expr_scaled <- scale(
  t(expr)
)

pca <- prcomp(
  expr_scaled,
  center = TRUE,
  scale. = TRUE
)

summary(pca)

plot(
  pca,
  type = "l",
  main = "Scree Plot"
)

# Latent components

clinical$Z1 <- pca$x[,1]
clinical$Z2 <- pca$x[,2]
clinical$Z3 <- pca$x[,3]

clinical$Z123 <- rowMeans(
  cbind(
    clinical$Z1,
    clinical$Z2,
    clinical$Z3
  )
)

# Kaplan-Meier groups

create_group <- function(x){
  
  ifelse(
    x > median(x, na.rm = TRUE),
    "High",
    "Low"
  )
  
}

clinical$Z1_group <- create_group(clinical$Z1)
clinical$Z2_group <- create_group(clinical$Z2)
clinical$Z3_group <- create_group(clinical$Z3)
clinical$Z123_group <- create_group(clinical$Z123)

# Kaplan-Meier analyses

fit_Z1 <- survfit(
  Surv(time_years,status) ~ Z1_group,
  data = clinical
)

fit_Z2 <- survfit(
  Surv(time_years,status) ~ Z2_group,
  data = clinical
)

fit_Z3 <- survfit(
  Surv(time_years,status) ~ Z3_group,
  data = clinical
)

fit_Z123 <- survfit(
  Surv(time_years,status) ~ Z123_group,
  data = clinical
)

# Unadjusted Cox models

cox_Z1 <- coxph(
  Surv(time_years,status) ~ Z1,
  data = clinical
)

cox_Z2 <- coxph(
  Surv(time_years,status) ~ Z2,
  data = clinical
)

cox_Z3 <- coxph(
  Surv(time_years,status) ~ Z3,
  data = clinical
)

cox_Z123 <- coxph(
  Surv(time_years,status) ~ Z123,
  data = clinical
)

# Multivariable Cox model

cox_model <- coxph(
  Surv(time_years,status) ~
    Z1 + Z2 + Z3 +
    age_z +
    subtype,
  data = clinical
)

summary(cox_model)

# Forest plot

coef_df <- data.frame(
  Variable = rownames(summary(cox_model)$coefficients),
  HR = exp(coef(cox_model)),
  Lower = exp(confint(cox_model)[,1]),
  Upper = exp(confint(cox_model)[,2])
)

# Proportional hazards assumption

ph_test <- cox.zph(
  cox_model
)

print(ph_test)

plot(ph_test)

# Time-dependent effects

cox_time <- coxph(
  Surv(time_years,status) ~
    Z123 +
    tt(Z123) +
    age_z +
    subtype,
  data = clinical,
  tt = function(x,t,...)
    x * log(t)
)

summary(cox_time)

# Nonlinearity assessment

cox_spline <- coxph(
  Surv(time_years,status) ~
    ns(Z123, df = 3) +
    age_z +
    subtype,
  data = clinical
)

anova(
  cox_model,
  cox_spline
)

# Tertile sensitivity analysis

clinical$Z123_tertile <- cut(
  clinical$Z123,
  breaks = quantile(
    clinical$Z123,
    probs = c(0,0.33,0.67,1),
    na.rm = TRUE
  ),
  labels = c(
    "Low",
    "Intermediate",
    "High"
  ),
  include.lowest = TRUE
)

fit_tertile <- survfit(
  Surv(time_years,status) ~ Z123_tertile,
  data = clinical
)

# Multiple testing correction

pvals <- summary(
  cox_model
)$coefficients[,5]

p.adjust(
  pvals,
  method = "BH"
)

# Subtype heterogeneity

summary(
  aov(
    Z1 ~ subtype,
    data = clinical
  )
)

summary(
  aov(
    Z2 ~ subtype,
    data = clinical
  )
)

summary(
  aov(
    Z3 ~ subtype,
    data = clinical
  )
)

kruskal.test(
  Z1 ~ subtype,
  data = clinical
)

kruskal.test(
  Z2 ~ subtype,
  data = clinical
)

kruskal.test(
  Z3 ~ subtype,
  data = clinical
)

# Subtype-specific Cox analyses

cox_luma <- coxph(
  Surv(time_years,status) ~
    Z123 + age_z,
  data = subset(
    clinical,
    subtype == "LumA"
  )
)

cox_lumb <- coxph(
  Surv(time_years,status) ~
    Z123 + age_z,
  data = subset(
    clinical,
    subtype == "LumB"
  )
)

cox_basal <- coxph(
  Surv(time_years,status) ~
    Z123 + age_z,
  data = subset(
    clinical,
    subtype == "Basal"
  )
)

cox_her2 <- coxph(
  Surv(time_years,status) ~
    Z123 + age_z,
  data = subset(
    clinical,
    subtype == "Her2"
  )
)

# Interaction model

cox_interaction <- coxph(
  Surv(time_years,status) ~
    Z123 * subtype +
    age_z,
  data = clinical
)

summary(cox_interaction)

# GO enrichment

loadings <- pca$rotation

get_top_genes <- function(pc){
  
  names(
    sort(
      loadings[,pc],
      decreasing = TRUE
    )[1:30]
  )
  
}

genes_Z1 <- get_top_genes(1)
genes_Z2 <- get_top_genes(2)
genes_Z3 <- get_top_genes(3)

ego_Z1 <- enrichGO(
  gene = genes_Z1,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)

ego_Z2 <- enrichGO(
  gene = genes_Z2,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)

ego_Z3 <- enrichGO(
  gene = genes_Z3,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)

dotplot(ego_Z1)
dotplot(ego_Z2)
dotplot(ego_Z3)

# E-value sensitivity analysis

Evalue <- function(HR){
  
  if(HR < 1)
    HR <- 1 / HR
  
  HR + sqrt(HR * (HR - 1))
}

EvalueCI <- function(LCL,UCL){
  
  if(LCL <= 1 & UCL >= 1)
    return(1)
  
  HR <- ifelse(
    LCL > 1,
    LCL,
    UCL
  )
  
  if(HR < 1)
    HR <- 1 / HR
  
  HR + sqrt(HR * (HR - 1))
}

Result <- data.frame(
  Variable = c("Z1","Z2","Z3"),
  HR = c(1.016,1.049,1.007),
  LCL = c(0.957,0.954,0.918),
  UCL = c(1.078,1.153,1.105)
)

Result$E_value <- round(
  sapply(
    Result$HR,
    Evalue
  ),
  3
)

Result$E_value_CI <- round(
  mapply(
    EvalueCI,
    Result$LCL,
    Result$UCL
  ),
  3
)

print(Result)

# Session information

sessionInfo()
