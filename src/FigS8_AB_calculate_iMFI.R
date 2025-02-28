

# setup -------------------------------------------------------------------

library(qs)
library(tidyverse)
library(data.table)
library(scattermore)
library(uwot)
library(ComplexHeatmap)
library(mixOmics)
select <- dplyr::select

hmcol <- circlize::colorRamp2(seq(0,1,length.out=8), rev(RColorBrewer::brewer.pal(8, "RdYlBu")))



# read data ---------------------------------------------------------------

# load integrated and annotated dataset
dt <- qs::qread("data/met_inhibitions/integrated_allLineages.qs") %>%
    dplyr::rename(TNF=TNFa)
dt$subset[dt$subset=="mDC_CD16+"] <- "mDC"


# FigS8_A -----------------------------------------------------------------


markers <- colnames(dt)[c(9:17, 19:26)]
msi <- dt |>
    summarise(across(all_of(markers), median), .by=subset) |>
    column_to_rownames(var="subset")

hm <- Heatmap(msi, name="MFI",
              col=hmcol, border=TRUE, rect_gp=gpar(col="black", lwd=.1),
              row_names_gp = gpar(fontsize=16/.pt),
              column_names_gp = gpar(fontsize=16/.pt),
              show_row_dend = FALSE, show_column_dend = FALSE)

cairo_pdf("outputs/FigS8_A_heatmap_subsets.pdf", width=4, height=2.5)
draw(hm)
dev.off()


# FigS8_B -----------------------------------------------------------------

set.seed(42)
dt.down <- dt %>%
    slice_sample(n=1e4, by=subset)
markers <- colnames(dt)[c(9:17, 19:26)]

# LDA initialization
lda.res <- lda(dt.down$subset ~ ., data=dt.down |> select(all_of(markers)))
coefs <- coef(lda.res)
lds <- as.matrix(dt.down |> select(all_of(markers))) %*% coefs

set.seed(42)
umap.res <- umap2(dt.down |> select(all_of(markers)),
                  fast_sgd=TRUE,
                  metric="cosine",
                  init=lds[,1:2],
                  n_threads=4, seed=42,
                  n_sgd_threads=4)

umap.res <- as.data.frame(umap.res)
umap.res$subset <- naturalsort::naturalfactor(dt.down$subset)
umap.res$lineage <- dt.down$lineage
umap.res$batch <- dt.down$batch
qs::qsave(umap.res, "data/FigS8_B_umap_coordinates.qs")

umap.res |>
    ggplot(aes(V1, V2, color=subset)) +
    geom_scattermore() +
    tidydr::theme_dr() +
    labs(x="LD-UMAP1", y="LD-UMAP2") +
    theme(panel.grid = element_blank()) +
    scale_color_manual(values = colorway::varibow(23))
ggsave("outputs/FigS8_B_umap_subsets.pdf", width=8, height=4, device=cairo_pdf)

# iMFI calculation --------------------------------------------------------

CUTOFF_IFNG = 0.15
CUTOFF_TNF  = 0.15

dt$IFNg_pos <- dt$IFNg >= CUTOFF_IFNG
dt$TNF_pos <- dt$TNF >= CUTOFF_TNF

# calculate counts of cytokine+ cells
count.cyto <-  dt |>
    summarise(
        IFNg = sum(IFNg_pos),
        TNF = sum(TNF_pos),
        .by=c(DonorID, Group, treatment, batch, lineage,  subset)
    )

# divide by counts of parent subsets
count.parent <- dt %>%
    count(DonorID, batch, treatment, subset, name="n_subset")

freq.parent <- left_join(count.parent, count.cyto,
                         by=c("DonorID", "batch", "treatment", "subset")) |>
    mutate(IFNg_freq = IFNg/n_subset*100,
           TNF_freq = TNF/n_subset*100) |>
    select(DonorID:TNF_freq)


# calculate MFI
mfi.parent <- dt |>
    summarise(
        IFNg_mfi = median(IFNg[IFNg_pos==TRUE]),
        TNF_mfi = median(TNF[TNF_pos==TRUE]),
        .by=c(DonorID, treatment, batch, lineage, subset)
    )


parent.freqmfi <- left_join(freq.parent, mfi.parent,
                            by=c("DonorID", "batch", "treatment", "subset", "lineage")) %>%
    mutate(IFNg_imfi = IFNg_freq * IFNg_mfi,
           TNF_imfi = TNF_freq * TNF_mfi) %>%
    mutate(Group = factor(Group, levels=c("RUR", "URB", "EUR"))) |>
    mutate(treatment = factor(treatment, levels=c("DMSO", "2DG", "6AN", "C75"))) %>%
    select(DonorID:subset, Group:lineage, IFNg_freq:TNF_imfi) %>%
    pivot_longer(c(IFNg_freq:TNF_imfi)) %>%
    mutate(type = str_split(name, "_", n=2, simplify=TRUE)[,2]) %>%
    mutate(cytokine = str_split(name, "_", n=2, simplify=TRUE)[,1]) %>%
    mutate(feature = paste(subset, cytokine, sep="_")) %>%
    select(DonorID:treatment, Group, lineage, subset, feature, type, value) %>%
    pivot_wider(names_from=type, values_from=value)


exclude.feature <- parent.freqmfi %>%
    summarise(mean = mean(freq), .by=feature) %>%
    filter(mean < 3.1) %>%
    pull(feature)
exclude.feature

parent.freqmfi <- parent.freqmfi %>% filter(!feature %in% exclude.feature)
qs::qsave(parent.freqmfi, "data/met_inhibitions/integrated_cytokine_allfeatures.qs")
