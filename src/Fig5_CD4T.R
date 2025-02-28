
suppressPackageStartupMessages({
    library(tidyverse)
})

PATHS <- list(
    EXV_limma = "outputs/EXV_limma_trend_test.rds",
    PMA_limma = "outputs/PMA_limma_trend_test.rds",
    MPL_limma = "outputs/MPL_limma_trend_test.rds"
)


# extract differential analysis results -----------------------------------

exv <- readRDS(PATHS$EXV_limma)
exv.coefs <- exv |>
    filter(lineage=="CD4T") |>
    rename(annotation=cluster) |>
    mutate(dataset = "EXV", .before="lineage") |>
    arrange(logFC) |>
    mutate(cluster = paste0("EXV.", c(2,1,6,13,5,14,3,15,4,12,20,16,18,11,10,19,9,17,8,7)), .before="annotation") |>
    select(dataset:lineage, cluster:annotation, linear_trend_coef=logFC, linear_trend_pvalue=pval_L, quadratic_trend_pvalue=fdr_Q, linear_trend_fdr=fdr_L) |>
    mutate(across(where(is.numeric), \(x) round(x, digits=5)))

exv.coefs |>
    filter(linear_trend_fdr  < .1 & quadratic_trend_pvalue > .1) |>
    mutate(cluster = fct_reorder(cluster, linear_trend_coef)) |>
    ggplot(aes(linear_trend_coef, cluster)) +
    geom_col()


pma <- readRDS(PATHS$PMA_limma)
pma.coefs <- pma |>
    filter(lineage=="CD4T") |>
    rename(annotation=cluster) |>
    mutate(dataset = "PMA", .before="lineage") |>
    arrange(logFC)|>
    filter(!annotation %in% c("CD4T_CD11c+CD161+CD27+TNF+", "CD4T_IL4IL5IL13+")) |>
    mutate(cluster = paste0("PMA.", c(18,17,7,19,5,15,3,16,14,6,1,12,10,11,8,4,2,13,9)),
           .before="annotation") |>
    mutate(cytokine_expressing = ifelse(str_detect(annotation, "IFN|TNF|IL"), "positive", "negative"), .after=annotation) |>
    select(dataset:annotation, linear_trend_coef=logFC,  linear_trend_pvalue=pval_L, quadratic_trend_pvalue=fdr_Q, linear_trend_fdr=fdr_L, cytokine_expressing) |>
    mutate(across(where(is.numeric), \(x) round(x, digits=6)))

pma.coefs |>
    filter(linear_trend_fdr  < .1 & quadratic_trend_pvalue > .1 & cytokine_expressing=="positive") |>
    mutate(cluster = fct_reorder(cluster, linear_trend_coef)) |>
    ggplot(aes(linear_trend_coef, cluster)) +
    geom_col()
