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

select <- dplyr::select

hmcol <- circlize::colorRamp2(seq(0,1,length.out=8), rev(RColorBrewer::brewer.pal(8, "RdYlBu")))
pals <- ggthemes::colorblind_pal()(8)
color_groups <- pals[c(2,3,4)]

conditions <- inhibitors <- c("2DG", "6AN", "C75")
groups <- c("RUR", "URB", "EUR")

source("src/stat_stars.R")
source("src/ordination_long.R")

theme_set(cowplot::theme_half_open(10))


# read data ---------------------------------------------------------------

parent.freqmfi <- qs::qread("data/met_inhibitions/integrated_cytokine_allfeatures.qs") |>
    filter(DonorID != "PK142")
ix <- which(parent.freqmfi$DonorID == "PK26" & parent.freqmfi$batch == "batch2")
parent.freqmfi <- parent.freqmfi[-ix, ]
markers <- unique(parent.freqmfi$feature)


# prepare data ------------------------------------------------------------

zscore <- parent.freqmfi %>%
    select(DonorID, treatment, batch, Group, feature, imfi=imfi) %>%
    pivot_wider(names_from=feature, values_from=imfi) |>
    mutate(treatment = factor(treatment, levels=c("DMSO", inhibitors)))

mat <- as.matrix(zscore[,markers])
mat <- mixOmics::impute.nipals(mat, 4)

modcombat <- model.matrix(~Group*treatment, data = zscore)
mat <- sva::ComBat(dat = t(mat),
                   batch = zscore$batch,
                   mean.only = FALSE,
                   mod = modcombat)
mat <- t(mat)
zscore[,markers] <- scale(mat)

groups <- c("RUR", "URB", "EUR")
results <- lapply(groups, \(g){
    zsub <- zscore |> filter(Group==g)

    sapply(inhibitors, \(cond){
        zsub2 <- zsub |> filter(treatment %in% c("DMSO", cond))

        rda.res <- rda(zsub2[,markers] ~ treatment, data=zsub2)
        pp <- scores(rda.res, scaling=1, choices=1:3, display="sites")
        set.seed(33); pvals <- anova.cca(rda.res, by="margin",
                                         permutations=how(nperm=5e3))
        pvals$`Pr(>F)`[1]
    })
})
names(results) <- groups

# stats -------------------------------------------------------------------

bind_rows(results, .id="group") |>
    mutate(group = factor(group, levels=groups)) |>
    pivot_longer(-1) |>
    mutate(name = factor(name, levels=conditions, labels=paste0(conditions, " vs. DMSO"))) |>
    ggplot(aes(group, -log10(value))) +
    geom_col(aes(fill=name), show.legend=FALSE) +
    geom_hline(yintercept=-log10(c(0.05, 0.01, 0.001)), lty=2) +
    geom_text(aes(label=paste0("P=",signif(value, 2))), vjust=1.5, cex=6/.pt) +
    cowplot::theme_half_open(10) +
    cowplot::panel_border(color="black") +
    labs(x=NULL, y="-Log10(P-value)") +
    scale_y_continuous(expand=expansion(mult=c(0, 0.1))) +
    facet_wrap2(~name, axes="all") +
    scale_fill_manual(values=RColorBrewer::brewer.pal(4, "Dark2")[-1]) +
    theme(strip.background = element_blank())
ggsave("outputs/Fig8_C_RDA_inhibitors_pvalues_byTreatment.pdf", width=4.5, height=2)
