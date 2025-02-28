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
    filter(str_detect(subset, "Monocytes")) |>
    rename(annotation=cluster) |>
    mutate(dataset = "EXV", .before="lineage") |>
    arrange(logFC)|>
    mutate(cluster = paste0("EXV.", c(1,2,3,8,6,5,7,4)), .before="annotation") |>
    select(dataset:lineage, cluster:annotation, linear_trend_coef=logFC, linear_trend_fdr=fdr_L, quadratic_trend_fdr=fdr_Q) |>
    mutate(across(where(is.numeric), \(x) round(x, digits=6)))

exv.coefs |>
    filter(linear_trend_fdr  < .1 & quadratic_trend_fdr > .1) |>
    mutate(cluster = fct_reorder(cluster, linear_trend_coef)) |>
    ggplot(aes(linear_trend_coef, cluster)) +
    geom_col()


pma <- readRDS(PATHS$PMA_limma)
pma.coefs <- pma |>
    filter(str_detect(cluster, "Monocytes")) |>
    rename(annotation=cluster) |>
    mutate(dataset = "PMA", .before="lineage") |>
    arrange(logFC)|>
    mutate(cluster = paste0("PMA.", c(8,4,2,7,6,5,1,3)), .before="annotation") |>
    mutate(cytokine_expressing = ifelse(str_detect(annotation, "IFN|TNF|IL"), "positive", "negative"), .after=annotation) |>
    select(dataset:lineage, cluster:annotation, linear_trend_coef=logFC, linear_trend_fdr=fdr_L, quadratic_trend_fdr=fdr_Q, cytokine_expressing) |>
    mutate(across(where(is.numeric), \(x) round(x, digits=6)))

pma.coefs |>
    filter(linear_trend_fdr  < .1 & cytokine_expressing=="positive") |>
    mutate(cluster = fct_reorder(cluster, linear_trend_coef)) |>
    ggplot(aes(linear_trend_coef, cluster)) +
    geom_col()



mpl <- readRDS(PATHS$MPL_limma)
mpl.coefs <- mpl |>
    filter(str_detect(cluster, "Monocytes")) |>
    rename(annotation=cluster) |>
    mutate(dataset = "MPL", .before="lineage") |>
    arrange(logFC) |>
    mutate(cluster = paste0("MPL", c(1,5,4,2,3)), .before="annotation") |>
    mutate(cytokine_expressing = ifelse(str_detect(annotation, "IFN|TNF|IL"), "positive", "negative"), .after=annotation) |>
    select(dataset:lineage, cluster:annotation, linear_trend_coef=logFC, linear_trend_fdr=fdr_L, quadratic_trend_fdr=fdr_Q, cytokine_expressing) |>
    mutate(across(where(is.numeric), \(x) round(x, digits=6)))

mpl.coefs |>
    filter(linear_trend_fdr  < .1 & quadratic_trend_fdr > .1 & cytokine_expressing=="positive") |>
    mutate(cluster = fct_reorder(cluster, linear_trend_coef)) |>
    ggplot(aes(linear_trend_coef, cluster)) +
    geom_col()
