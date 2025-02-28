# functions --------------------
suppressPackageStartupMessages({
    library(ggforce)
    library(limma)
    library(tidyverse)
    library(qs)
    library(scattermore)
    library(ggbeeswarm)
    library(cowplot)
    library(ggpubr)
    library(data.table)
    library(patchwork)
    library(ComplexHeatmap)


    source("src/utils.R")
    source("src/utils_plot.R")
    source("src/weightedcontrasts.R")

})


PATHS <- list()
PATHS$data <- "data/EXV_sce_annotated.qs"
PATHS$counts <- "data/EXV_counts_subclusters.qs"
PATHS$results <- "outputs/"
PATHS$metadata <- "data/PMA_metadata.rds"

# prepare data -------------------------------------

sce <- qs::qread(PATHS$data, nthreads=6)
markersALL <- rownames(sce)
sce.dt <- sce2dt(sce)

# save cluster medians
cluster.medians <- sce.dt[, lapply(.SD, median), .SDcols = markersALL, by = cluster]
qs::qsave(cluster.medians, "data/EXV_cluster_median_markers.qs")


## calculate frequencies ---------------------------------------

countsALL <- qs::qread(PATHS$counts)
metaLineage <- countsALL %>%
    select(lineage, subset, cluster) %>%
    distinct()

mat <- countsALL %>%
    mutate(freq = asin(sqrt(n_cluster / n_lineage))) %>%
    select(patient_id:sex, cluster, freq) %>%
    pivot_wider(names_from = cluster, values_from = freq, values_fill = 0) %>%
    mutate(group = factor(group2, levels = c("RUR", "URB", "EUR"), ordered = F), .after = group2) %>%
    # remove reference samples
    filter(!grepl("REF", patient_id)) %>%
    filter(!grepl("LD413", patient_id))

freqLineage <- countsALL %>%
    # remove reference samples
    filter(!grepl("REF", patient_id)) %>%
    filter(!grepl("LD413", patient_id)) %>%
    select(patient_id, group2, barcode_id, age, sex, lineage, n_total, n_lineage) %>%
    distinct() %>%
    mutate(freqLineage = n_lineage/n_total*100)


matMeta <- mat %>% select(patient_id:sex)
matCounts <- t(mat[, 7:ncol(mat)])
colnames(matCounts) <- matMeta$patient_id
clusterMeta <- distinct(as.data.table(colData(sce)[, c("lineage", "subset", "cluster")]))


# trend test limma --------------------------------------------------------

matMeta <- mat %>% dplyr::select(patient_id:sex)
contrasts(matMeta$group2) <- contr.poly.weighted(matMeta$group2)
design <- model.matrix(~ group2 + age + sex + barcode_id, data = matMeta)
fit <- lmFit(matCounts, design, method = "robust", maxit = 100)
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
    left_join(clusterMeta, ., by = "cluster")
res.LQ.weighted


out <- paste0(PATHS$results, "EXV_limma_trend_test.rds")
saveRDS(res.LQ.weighted, out)
