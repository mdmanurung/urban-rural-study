# setup -------------------------------------------------------------------

library(qs)
library(tidyverse)
library(data.table)
library(scattermore)
library(uwot)
library(ComplexHeatmap)
library(mixOmics)
library(patchwork)
library(vegan)
library(ggh4x)
library(ggbeeswarm)
library(cowplot)
library(ConsensusClusterPlus)

select <- dplyr::select

hmcol <- circlize::colorRamp2(seq(0,1,length.out=8), rev(RColorBrewer::brewer.pal(8, "RdYlBu")))
pals <- ggthemes::colorblind_pal()(8)
color_groups <- pals[c(2,3,4)]
inhibitors <- c("C75", "2DG", "6AN")

source("src/stat_stars.R")

theme_set(theme_half_open(10))


# read data ---------------------------------------------------------------

parent.freqmfi <- qs::qread("data/met_inhibitions/integrated_cytokine_allfeatures.qs") |>
    filter(DonorID != "PK142")

ix <- which(parent.freqmfi$DonorID == "PK26" & parent.freqmfi$batch == "batch2")
parent.freqmfi <- parent.freqmfi[-ix, ]



# FigS8_C -----------------------------------------------------------------

zscore <- parent.freqmfi %>%
    select(DonorID, treatment, feature, imfi) %>%
    mutate(imfi = as.numeric(scale(imfi)), .by=c(DonorID, feature)) %>%
    pivot_wider(names_from=feature, values_from=imfi)


mat <- as.matrix(t(zscore[,3:ncol(zscore)]))
mat[is.na(mat)] <- 0
cons <- ConsensusClusterPlus(t(mat), maxK=10, reps=100, pItem=.8,
                             innerLinkage="ward.D2", finalLinkage="ward.D2", distance="spearman")
clust <- cons[[5]]$consensusClass

hm <- Heatmap(mat,
              name="Z-score",
              column_split=factor(zscore$treatment, levels=c("DMSO", "C75", "2DG", "6AN")),
              cluster_column_slices = FALSE,
              row_split = clust,
              clustering_method_rows = "ward.D2", clustering_distance_rows="euclidean",
              show_row_dend = FALSE, show_column_dend = FALSE,
              row_names_gp = gpar(fontsize=8),
              border=TRUE)

cairo_pdf("outputs/FigS8_C_heatmap_zscores.pdf", width=6, height=5)
draw(hm)
dev.off()



# Fig8_B ------------------------------------------------------------------
# PCA on z-scores

zscore <- parent.freqmfi %>%
    select(DonorID, treatment, batch, Group, feature, imfi) %>%
    pivot_wider(names_from=feature, values_from=imfi) |>
    mutate(treatment = factor(treatment, levels=c("DMSO", inhibitors)))


mat <- as.matrix(t(zscore[,5:ncol(zscore)]))
mat.corr <- limma::removeBatchEffect(mat,
                                     batch = zscore$batch,
                                     group = with(zscore, paste(treatment, Group)))
zscore.corr <- zscore
zscore.corr[,5:ncol(zscore.corr)] <- t(mat.corr)

pc.res <- mixOmics::pca(zscore.corr[,5:ncol(zscore.corr)], scale=TRUE, center=T, multilevel=zscore.corr$DonorID)
p1 <- plotIndiv(pc.res, group=as.character(zscore.corr$treatment), legend=TRUE, star=TRUE, ellipse=FALSE, title="PCA", pch=16,
                col.per.group = RColorBrewer::brewer.pal(4, "Dark2"))

plotIndiv(pc.res, group=zscore.corr$treatment, legend=TRUE, star=TRUE, ellipse=FALSE, title="PCA")


zscore.corr$PC1 <- pc.res$x[,1]
zscore.corr$PC2 <- pc.res$x[,2]
p1 <- zscore.corr |>
    ggplot(aes(PC1, PC2, color=treatment)) +
    geom_point(aes(shape=treatment), cex=1.5) +
    stat_stars() +
    stat_ellipse() +
    cowplot::theme_half_open(10) +
    scale_color_brewer(palette="Dark2") +
    theme(aspect.ratio = 22/38) +
    labs(x="PC1 (38% var. expl.)", y="PC2 (22% var. expl.)")
p1
ggsave("outputs/Fig8_B_pca_inhibitors_imfi_treatment.pdf", width=3, height=2, device=cairo_pdf, scale=1.5)


