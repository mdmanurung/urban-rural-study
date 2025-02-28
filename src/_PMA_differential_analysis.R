# functions --------------------
suppressPackageStartupMessages({
    library(tidyverse)
    library(qs)
    library(scattermore)
    library(ggbeeswarm)
    library(cowplot)
    library(ggpubr)
    library(data.table)
    library(ComplexHeatmap)
    library(patchwork)
    library(limma)
    library(ggforce)

    source("src/utils.R")
    source("src/utils_plot.R")
    source("src/weightedcontrasts.R")
})

PATHS <- list(
    data = "data/PMA_data.qs",
    meta = "data/PMA_metadata.rds",
    results = "outputs/"
)

# read data ---------------------------------------------------------------

dat <- qs::qread(PATHS$data, nthreads=6)
markersAll <- colnames(dat[,CD45:IL10])

cluster.medians <- dat[, lapply(.SD, median), .SDcols = markersAll, by = subclusterPG_merged]
qs::qsave(cluster.medians, "data/PMA_cluster_median_markers.qs")


# prepare count tables ----------------------------------------------------

countsLineage <- dplyr::count(dat, patient_id, lineage, name = "n_lineage")
countsCluster <- dplyr::count(dat, patient_id, subclusterPG_merged, name = "n_cluster")
metaLineage <- dat %>%
    select(lineage, subclusterPG_merged) %>%
    distinct()


metaSample <- readRDS(PATHS$meta)
counts <- left_join(metaLineage, countsCluster, by = "subclusterPG_merged") %>%
    left_join(countsLineage, ., by = c("patient_id", "lineage")) %>%
    left_join(metaSample, ., by = "patient_id") %>%
    as.data.table() %>%
    mutate(freq_cluster_perLineage = n_cluster / n_lineage * 100) %>%
    mutate(freq_cluster_perLineage_trans = asin(sqrt(freq_cluster_perLineage / 100))) %>%
    drop_na() %>%
    # remove reference sample
    filter(patient_id != "LD413")
qs::qsave(counts, "data/PMA_counts_subclusters.qs")

matCounts <- counts %>%
    select(patient_id, subclusterPG_merged, freq_cluster_perLineage_trans) %>%
    pivot_wider(names_from = patient_id, values_from = freq_cluster_perLineage_trans, values_fill = 0) %>%
    column_to_rownames(var = "subclusterPG_merged")

matMeta <- counts %>%
    select(patient_id, group3, age, sex, batch) %>%
    mutate(group3 = factor(group3, levels=c("RUR", "URB", "EUR"))) |>
    distinct()


# differential analysis limma ---------------------------------------------

contrasts(matMeta$group3) <- contr.poly.weighted(matMeta$group3)
designTrend <- model.matrix(~ group3 + age + sex + batch, data = matMeta)
fit <- lmFit(matCounts, designTrend, method = "robust", maxit = 100)
efit <- eBayes(fit, robust = TRUE)
res.L <- topTable(efit, coef = 2, number = Inf) %>%
    as.data.frame() %>%
    rownames_to_column(var = "cluster") %>%
    as_tibble()%>%
    dplyr::rename(pval_L = P.Value, fdr_L = adj.P.Val)
res.Q <- topTable(efit, coef = 3, number = Inf) %>%
    as.data.frame() %>%
    rownames_to_column(var = "cluster") %>%
    as_tibble() %>%
    dplyr::rename(pval_Q = P.Value, fdr_Q = adj.P.Val)
res.LQ.weighted <- left_join(res.L,
                             res.Q %>% select(cluster, pval_Q, fdr_Q),
                             by = "cluster") %>%
    left_join(metaLineage %>% dplyr::rename(cluster=2), ., by = "cluster")

out <- paste0(PATHS$results, "PMA_limma_trend_test.rds")
saveRDS(res.LQ.weighted, out)
