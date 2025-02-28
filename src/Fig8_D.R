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
library(emmeans)
require(nlme)

select <- dplyr::select

hmcol <- circlize::colorRamp2(seq(0,1,length.out=8), rev(RColorBrewer::brewer.pal(8, "RdYlBu")))
pals <- ggthemes::colorblind_pal()(8)
color_groups <- pals[c(2,3,4)]
inhibitors <- c("2DG", "6AN", "C75")

source("src/stat_stars.R")

theme_set(theme_half_open(10))
conditions <- c("2DG", "6AN", "C75")


# read data ---------------------------------------------------------------

parent.freqmfi <- qs::qread("data/met_inhibitions/integrated_cytokine_allfeatures.qs") |>
    filter(DonorID != "PK142")

ix <- which(parent.freqmfi$DonorID == "PK26" & parent.freqmfi$batch == "batch2")
parent.freqmfi <- parent.freqmfi[-ix, ]


# stats -------------------------------------------------------------------

zscore <- parent.freqmfi %>%
    select(DonorID, treatment, batch, Group, feature, imfi) %>%
    pivot_wider(names_from=feature, values_from=imfi) |>
    mutate(treatment = factor(treatment, levels=c("DMSO", inhibitors)))

zscore.nest <- zscore |>
    pivot_longer(where(is.numeric)) |>
    group_by(name) |>
    nest()

features <- zscore.nest$name

fits <- list()
emms <- list()
for (i in 1:nrow(zscore.nest)){
    y <-zscore.nest$data[[i]] %>% drop_na(value)

    fit <- lme(scale(value) ~ Group*treatment,
               data=y,
               random=~1|batch/DonorID,
               weights=varIdent(form=~1|Group*treatment),
               method="REML",
               control = lmeControl(opt = "optim"))

    fits[[i]] <- fit
    emms[[i]] <- data.frame(emmeans(fit, trt.vs.ctrl ~ treatment|Group, data=y, adjust="none")$contrasts)
}

names(emms) <- features
emms_all <- bind_rows(emms, .id="subset") %>%
    mutate(contrast = str_remove_all(contrast, " - DMSO")) %>%
    mutate(fdr = p.adjust(p.value), .by=contrast)



# Fig8_D: regression_coef ---------------------------------------------------------


FDR_THRESHOLD <- 0.1

sig.list <- list()
plots <- list()
for (inhib in inhibitors){
    emms_sub <- emms_all %>% filter(contrast==inhib)

    list.sig <- lapply(c("RUR", "URB", "EUR"), function(x){
        emms_sub %>%
            filter(Group == x) %>%
            filter(fdr < FDR_THRESHOLD) %>%
            pull(subset)
    })
    sigs <- unique(Reduce(c, list.sig))
    sig.list[[inhib]] <- sigs

    ords <-emms_sub %>%
        filter(subset %in% sigs) %>%
        filter(Group == "EUR") %>%
        mutate(subset = fct_reorder(subset, estimate))

    plots[[inhib]] <- emms_sub %>%
        filter(subset %in% sigs) %>%
        mutate(direction = case_when(
            fdr < FDR_THRESHOLD & estimate >0 ~ "Up",
            fdr < FDR_THRESHOLD & estimate <0 ~ "Down",
            TRUE ~ "NS"
        )) |>
        mutate(subset = factor(subset, levels=levels(ords$subset))) %>%
        ggplot(aes(y=subset, x=estimate, fill=direction)) +
        ggstats::geom_stripped_rows(odd="white", even="gray90") +
        geom_vline(xintercept=0, lty=2) +
        geom_col() +
        labs(title=paste0(inhib, " vs. DMSO")) +
        scale_fill_manual(values=c("Up" = "firebrick", "NS" = "grey60", "Down" = "steelblue2")) +
        cowplot::theme_cowplot(8) +
        theme(axis.line.y = element_blank()) +
        scale_y_discrete(position="right") +
        labs(x = "Regression Coef.", y=NULL) +
        facet_wrap2(~Group,
                    strip = strip_themed(
                        background_x = list(
                            element_rect(fill = color_groups),
                            element_rect(fill = color_groups[2]),
                            element_rect(fill = color_groups[3])
                        )
                    ))
}
wrap_plots(plots, nrow=3, heights = sapply(sig.list, length)) + plot_layout(guides="collect")
ggsave("outputs/Fig8_D_inhibitor_coefs.pdf", height=8, width=5)



# FigS8_D -----------------------------------------------------------------

p1 <- zscore |>
    mutate(NK_CD56dim_IFNg = log2(NK_CD56dim_IFNg / NK_CD56dim_IFNg[treatment=="DMSO"]), .by=DonorID) |>
    filter(treatment != "DMSO") |>
    mutate(treatment = factor(treatment, levels=conditions)) |>

    ggplot(aes(Group, NK_CD56dim_IFNg)) +
    geom_hline(yintercept=0, lty=2) +
    geom_boxplot(outlier.shape=NA) +
    geom_quasirandom(aes(color=Group)) +
    facet_wrap(~treatment) +
    labs(x=NULL, y="LFC vs. DMSO", title="IFNg+ CD56dim NK") +
    scale_color_manual(values=color_groups) +
    geom_text(
        data = emms_all |> filter(subset=="NK_CD56dim_IFNg") |> mutate(treatment=factor(contrast, levels=conditions)),
        inherit.aes=FALSE,
        aes(x=Group, y=Inf, label=signif(fdr, 2)),
        vjust=1, cex=8/.pt
    )


p2 <- zscore |>
    mutate(CD4T_EM_IFNg = log2(CD4T_EM_IFNg / CD4T_EM_IFNg[treatment=="DMSO"]), .by=DonorID) |>
    filter(treatment != "DMSO") |>
    mutate(treatment = factor(treatment, levels=conditions)) |>

    ggplot(aes(Group, CD4T_EM_IFNg)) +
    geom_hline(yintercept=0, lty=2) +
    geom_boxplot(outlier.shape=NA) +
    geom_quasirandom(aes(color=Group)) +
    facet_wrap(~treatment) +
    labs(x=NULL, y="LFC vs. DMSO", title="IFNg+ EM CD4+T") +
    scale_color_manual(values=color_groups) +
    geom_text(
        data = emms_all |> filter(subset=="CD4T_EM_IFNg") |> mutate(treatment=factor(contrast, levels=conditions)),
        inherit.aes=FALSE,
        aes(x=Group, y=Inf, label=signif(fdr, 2)),
        vjust=1, cex=8/.pt
    )


p3 <- zscore |>
    mutate(Monocytes_Classical_TNF = log2(Monocytes_Classical_TNF / Monocytes_Classical_TNF[treatment=="DMSO"]), .by=DonorID) |>
    filter(treatment != "DMSO") |>
    mutate(treatment = factor(treatment, levels=conditions)) |>

    ggplot(aes(Group, Monocytes_Classical_TNF)) +
    geom_hline(yintercept=0, lty=2) +
    geom_boxplot(outlier.shape=NA) +
    geom_quasirandom(aes(color=Group)) +
    facet_wrap(~treatment) +
    labs(x=NULL, y="LFC vs. DMSO", title="TNF+ Classical Monocytes") +
    scale_color_manual(values=color_groups) +
    geom_text(
        data = emms_all |> filter(subset=="Monocytes_Classical_TNF") |> mutate(treatment=factor(contrast, levels=conditions)),
        inherit.aes=FALSE,
        aes(x=Group, y=Inf, label=signif(fdr, 2)),
        vjust=1, cex=8/.pt
    )

p1/p2/p3 + plot_layout(guides="collect")
ggsave("outputs/FigS8_D.pdf", width=6, height=6)
