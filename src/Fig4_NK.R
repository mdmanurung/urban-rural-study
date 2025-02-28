
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
    rename(annotation=cluster) |>
    filter(str_detect(annotation, "NK")) |>
    mutate(dataset = "EXV", .before="lineage") |>
    arrange(annotation) |>
    arrange(logFC) |>
    mutate(cluster = paste0("EXV.", c(12,3,15,1,10,11,5,7,13,2,9,14,6,8,4)), .before="annotation") |>
    select(dataset:lineage, cluster:annotation, linear_trend_coef=logFC,  linear_trend_fdr=fdr_L, quadratic_trend_fdr=fdr_Q) |>
    mutate(across(where(is.numeric), \(x) round(x, digits=4)))

exv.coefs |>
    filter(linear_trend_fdr  < .1 & quadratic_trend_fdr > .1) |>
    mutate(cluster = fct_reorder(cluster, linear_trend_coef)) |>
    ggplot(aes(linear_trend_coef, cluster)) +
    geom_col()


pma <- readRDS(PATHS$PMA_limma)
pma.coefs <- pma |>
    filter(lineage %in% c("NK")) |>
    rename(annotation=cluster) |>
    mutate(dataset = "PMA", .before="lineage") |>
    mutate(cluster = paste0("PMA.", row_number()), .before="annotation") |>
    # arrange(annotation) |>
    arrange(logFC)|>
    mutate(cytokine_expressing = ifelse(str_detect(annotation, "IFN|TNF|IL"), "positive", "negative"), .after=annotation) |>
    select(dataset:lineage, cluster:annotation, linear_trend_coef=logFC, linear_trend_fdr=fdr_L, quadratic_trend_pvalue=fdr_Q, cytokine_expressing) |>
    mutate(across(where(is.numeric), \(x) round(x, digits=4)))

pma.coefs |>
    filter(linear_trend_fdr < .1 & quadratic_trend_pvalue > .01 & cytokine_expressing == "positive") |>
    ggplot(aes(linear_trend_coef, cluster)) +
    geom_col()


mpl <- readRDS(PATHS$MPL_limma)
mpl.coefs <- mpl |>
    filter(lineage %in% c("NK")) |>
    rename(annotation=cluster) |>
    mutate(dataset = "MPL", .before="lineage") |>
    arrange(annotation) |>
    arrange(logFC)|>
    mutate(cluster = paste0("MPL.", c(13,11,14,15,17,12,2,10,9,4,5,8,16,1,3,7,6)), .before="annotation") |>
    mutate(cytokine_expressing = ifelse(str_detect(annotation, "IFN|TNF|IL"), "positive", "negative"), .after=annotation) |>
    select(dataset:lineage, cluster:annotation, linear_trend_coef=logFC, linear_trend_fdr=fdr_L, quadratic_trend_fdr=fdr_Q, cytokine_expressing) |>
    mutate(across(where(is.numeric), \(x) round(x, digits=6)))

mpl.coefs |>
    filter(linear_trend_fdr < .1 & cytokine_expressing == "positive") |>
    ggplot(aes(linear_trend_coef, cluster)) +
    geom_col()
