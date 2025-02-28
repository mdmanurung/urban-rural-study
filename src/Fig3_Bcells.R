

suppressPackageStartupMessages({
    library(tidyverse)
})

PATHS <- list(
    EXV_limma = "outputs/EXV_limma_trend_test.rds",
    PMA_limma = "outputs/PMA_limma_trend_test.rds"
)


# extract differential analysis results -----------------------------------

exv <- readRDS(PATHS$EXV_limma)

# effect size
exv.coefs <- exv |>
    filter(lineage=="B") |>
    rename(annotation=cluster) |>
    mutate(dataset = "EXV", .before="lineage") |>
    mutate(cluster = paste0("EXV.", row_number()), .before="annotation") |>
    arrange(logFC) |>
    select(dataset:annotation, trend_test_coef=logFC, pvalue_linear=pval_L, pvalue_quadratic=pval_Q, fdr_linear=fdr_L) |>
    mutate(across(where(is.numeric), \(x) round(x, digits=4))) |>
    mutate(cluster = fct_reorder(cluster, trend_test_coef))
exv.coefs

exv.coefs |>
    filter(fdr_linear < .1 & pvalue_quadratic > .1) |>
    ggplot(aes(trend_test_coef, cluster)) +
    geom_col()


pma <- readRDS(PATHS$PMA_limma)
pma.coefs <- pma |>
    filter(lineage=="B") |>
    rename(annotation=cluster) |>
    mutate(dataset = "PMA", .before="lineage") |>
    mutate(cluster = paste0("PMA.", row_number()), .before="annotation") |>
    arrange(logFC)|>
    mutate(cytokine_expressing = ifelse(str_detect(annotation, "IFN|TNF|IL"), "positive", "negative"), .after=annotation) |>
    select(dataset:annotation, trend_test_coef=logFC,  pvalue_linear=pval_L, pvalue_quadratic=pval_Q, fdr_linear=fdr_L, cytokine_expressing) |>
    mutate(across(where(is.numeric), \(x) round(x, digits=4)))

pma.coefs |>
    filter(fdr_linear < .1 & pvalue_quadratic > .1 & cytokine_expressing == "positive") |>
    ggplot(aes(trend_test_coef, cluster)) +
    geom_col()
