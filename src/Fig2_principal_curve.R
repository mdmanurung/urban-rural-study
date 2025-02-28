suppressPackageStartupMessages({
    library(tidyverse)
    library(here)
    library(glue)
    library(patchwork)
    library(Hmisc)
    library(nlme)
    library(emmeans)
    library(PMCMRplus)
    library(princurve)
    library(ggpubr)
    library(ggbeeswarm)
    library(vegan)
    library(ggnewscale)
})

PATHS <- list(
    latent_factors = "outputs/mofa/final_mofa_model_factors.rds",
    meta = "data/PMA_metadata.rds",
    mofa_model = "outputs/mofa/final_mofa_model.rds"
)

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflicts_prefer(stats::predict)

color_groups <- c("#E69F00", "#56B4E9", "#009E73")
color_groups <- purrr::set_names(color_groups, c("RUR", "URB", "EUR"))
theme_set(theme_pubr(10, legend="right"))



# load data ---------------------------------------------------------------

meta <- readRDS(PATHS$meta) %>%
    rename(sample = patient_id) %>%
    as.data.frame() %>%
    mutate(group3 = factor(group3, levels=c("RUR", "URB", "EUR")))
rownames(meta) <- meta$sample

# load mofa latent factors
latent.factors <- readRDS(PATHS$latent_factors)

# plot latent factors -----------------------------------------------------


p <- lapply(c(1:8), function(i){

    message(i)
    fct <- paste0("Factor", i)

    # pairwise p-values
    pairs <- PMCMRplus::kwAllPairsConoverTest(x=latent.factors[[fct]], latent.factors$group3)
    pairs.df <- data.frame(xmin=c("RUR", "RUR", "URB"),
                           xmax=c("URB", "EUR", "EUR"),
                           p.value=scales::pvalue(c(pairs$p.value[1,1], pairs$p.value[2,1], pairs$p.value[2,2]), accuracy=1e-4),
                           y.position=max(latent.factors[[fct]])*1.1)

    # trend p values
    trend <- PMCMRplus::jonckheereTest(latent.factors[[fct]], latent.factors$group3)
    trend.pval <- scales::pvalue(trend$p.value, accuracy=1e-4)

    # errorbars
    fit <- gls(as.formula(paste0(fct, " ~ group3")), data=latent.factors, weights=varIdent(form=~1|group3))
    means <- data.frame(emmeans(fit, pairwise ~ group3, mode="df.error")$emmeans)

    # plots
    latent.factors |>
        ggplot(aes(group3, .data[[fct]])) +
        geom_quasirandom(aes(color=group3)) +
        geom_errorbar(data=means, aes(x=group3, ymin=lower.CL, ymax=upper.CL), inherit.aes=FALSE, width=.5) +
        geom_point(data=means, aes(x=group3, y=emmean, shape=group3, fill=group3), inherit.aes=FALSE, cex=3) +
        ggpubr::geom_bracket(data=pairs.df, step.increase=.08, tip.length=0.01, label.size=7/.pt) +
        labs(x=NULL, y=fct, subtitle=glue::glue("Trend P={trend.pval}")) +
        # coord_cartesian(clip="off") +
        scale_fill_manual(values=color_groups) +
        scale_color_manual(values=color_groups) +
        scale_shape_manual(values=c(21,22,23)) +
        coord_cartesian(clip="off") +
        theme(legend.position = "none")
})
p <- wrap_plots(p, nrow=2)
ggsave("outputs/FigS2_A.pdf", p, dpi=600, width=8, height=4.5)



# Multidimensional scaling and principal curve ---------------------------------------------------------------------


latent.factors <- readRDS("outputs/mofa/final_mofa_model_factors.rds")
latent.factors.coords <- latent.factors[, grepl("Factor", colnames(latent.factors))]

mds.res <- cmdscale(dist(latent.factors.coords))
latent.factors$mds1 <- mds.res[,1]
latent.factors$mds2 <- mds.res[,2]

set.seed(42)
pcurve <- principal_curve(as.matrix(latent.factors[,c("mds1", "mds2")]), maxit=100, thresh = 0.0005)
latent.factors$lambda <- pcurve$lambda

lf.reord <- latent.factors
ix <- order(lf.reord$lambda)
pcurve.coords <- pcurve$s[ix,]
lambda <- lf.reord$lambda[ix]


adonis <- adonis2(lf.reord[,c("mds1", "mds2")] ~ group3, data=lf.reord, method="euclidean", permutations=1e4)
adonis # PERMANOVA pvalues

p1 <- lf.reord |>
    ggplot(aes(mds1, mds2, color=group3)) +
    geom_point(cex=3, aes(shape=group3, color=group3, fill=group3)) +
    stat_stars(show.legend = FALSE) +
    stat_conf_ellipse(show.legend = FALSE) +
    tidydr::theme_dr() +
    labs(x="MDS1", y="MDS2") +
    scale_fill_manual(values=color_groups) +
    scale_color_manual(values=color_groups) +
    scale_shape_manual(values=c(21,22,23)) +
    theme(panel.grid = element_blank(),
          legend.position = c("inside"),
          legend.position.inside = c(1,1),
          legend.justification = c(1,1)) +
    coord_equal()

p2 <- lf.reord |>
    ggplot(aes(mds1, mds2, color=group3)) +
    geom_point(cex=3, show.legend=FALSE) +
    geom_segment(x=lf.reord$mds1[ix], xend=pcurve.coords[,1],
                 y=lf.reord$mds2[ix], yend=pcurve.coords[,2],
                 inherit.aes=FALSE, color="gray80", show.legend=FALSE) +
    scale_fill_manual(values=color_groups) +
    scale_color_manual(values=color_groups) +
    scale_shape_manual(values=c(21,22,23)) +

    new_scale_color() +
    geom_path(data=pcurve.coords, aes(mds1, mds2, color=lambda), linewidth=2, inherit.aes=FALSE) +
    scale_color_viridis_c(name="Lambda", limits=c(0, 11), breaks=c(0,11), labels=c("min", "max")) +
    tidydr::theme_dr() +
    theme(panel.grid = element_blank(),
          legend.position = c("inside"),
          legend.position.inside = c(1,1),
          legend.justification = c(1,1)) +
    labs(x="MDS1", y="MDS2") +
    guides(color=guide_colorbar(direction="horizontal", legend.ticks=element_blank(), title.position="top")) +

    coord_equal()

fit <- gls(rank(lambda) ~ group3, data=latent.factors)
emm.trend <- emmeans(fit, poly ~ group3, adjust="none")
(trend.pval <- as.data.frame(emm.trend$contrasts)[1,"p.value"]) # lambda trend-test P-value

emm.pairs <- emmeans(fit, pairwise ~ group3)
pvals <- as.data.frame(emm.pairs$contrasts) |>
    separate(contrast, " - ", into=c("xmin", "xmax")) |>
    mutate(y.position = max(latent.factors$lambda)) |>
    mutate(p.value = rstatix::p_round(p.value))

fit <- gls(lambda ~ group3, data=latent.factors)
means <- data.frame(emmeans(fit, pairwise ~ group3)$emmeans)

p3 <- latent.factors |>
    ggplot(aes(group3, lambda)) +
    geom_quasirandom(aes(color=group3)) +
    geom_errorbar(data=means, aes(x=group3, ymin=lower.CL, ymax=upper.CL), inherit.aes=FALSE, width=.5) +
    geom_point(data=means, aes(x=group3, y=emmean, shape=group3, fill=group3), inherit.aes=FALSE, cex=3) +
    ggpubr::geom_bracket(data=pvals, step.increase=.08, tip.length=0.01) +
    labs(x=NULL, y="Lambda", subtitle="Trend P < 0.0001") +
    # coord_cartesian(clip="off") +
    scale_fill_manual(values=color_groups) +
    scale_color_manual(values=color_groups) +
    scale_shape_manual(values=c(21,22,23)) +
    theme(legend.position = "none")



p1+p2
ggsave("outputs/Fig2_CD_MDSandPrinCurve.pdf", dpi=600, width=8, height=3)

p3
ggsave("outputs/Fig2_D_plotsLambda.pdf", dpi=600, width=2, height=2.5)


# Correlate features with lambda ------------------------------------------

mofa <- readRDS(PATHS$mofa_model)
mofa.data <- get_data(mofa, as.data.frame=TRUE)
mofa.df <- mofa.data |>
    select(feature:value) |>
    pivot_wider(names_from=feature, values_from=value) |>
    mutate(sample = as.character(sample)) |>
    column_to_rownames(var="sample")

results <- list()
for (i in 1:ncol(mofa.df)){
    cors <- cor.test(mofa.df[,i], latent.factors$lambda, method="spearman", continuity=TRUE, exact=FALSE)
    results[[i]] <- as.numeric(c(cors$p.value, cors$estimate))
}
df <- data.frame(do.call(rbind, results))
names(df) <- c("p.value", "rho")
df$feature <- colnames(mofa.df)
df$padj <- qvalue::qvalue(df$p.value)$qvalues


df.cors <- df |>
    left_join(mofa.data |> select(view, feature) |> distinct(), by="feature") |>
    mutate(
        lineage = case_when(
            str_detect(feature, "CD4T|Th2|Th17|Treg|Th1") ~ "CD4T",
            str_detect(feature, "^B|Bcells|IgG") ~ "B cells",
            str_detect(feature, "CD8T|Tc2") ~ "CD8T",
            str_detect(feature, "gdT") ~ "gdT",
            str_detect(feature, "DPT|DNT|NKT") ~ "Unconv. T",
            str_detect(feature, "NK|ILC") ~ "NK/ILC",
            str_detect(feature, "Monocytes|mDC|pDC") ~ "Myeloids",
            TRUE ~ "Misc"
        )
    )

selected <- df.cors |>
    filter(view != "metflow") |>
    mutate(direction = sign(rho)) |>
    filter(padj < .1) |>
    filter(feature != "mDC") |>
    filter(lineage %in% c("B cells", "NK/ILC", "CD4T", "Myeloids")) |>
    mutate(lineage = factor(lineage, levels=c("B cells", "NK/ILC", "CD4T", "Myeloids"))) |>
    slice_max(order_by = rho, n=1, by=c(direction, lineage), with_ties=FALSE)

df.cors2 <- df.cors |>
    filter(lineage %in% c("B cells", "NK/ILC", "CD4T", "Myeloids")) |>
    mutate(lineage = factor(lineage, levels=c("B cells", "NK/ILC", "CD4T", "Myeloids"))) |>
    filter(view != "metflow") |>
    mutate(direction = factor(sign(rho)))

p <- df.cors2 |>
    ggplot(aes(rho, -log10(padj))) +
    geom_point() +
    geom_point(data = df.cors2 |>  filter(padj < .1), aes(color=direction))  +
    facet_wrap(~lineage, nrow=1) +
    ggrepel::geom_text_repel(data=selected, aes(label=feature), size=7/.pt) +
    labs(x="Spearman's Rho(Lambda, Features)", y="-log10(q-value)") +
    geom_hline(yintercept=1, lty=2) +
    cowplot::theme_cowplot(10) +
    cowplot::panel_border(color="black") +
    theme(strip.background = element_blank(), strip.text = element_text(face="bold"),
          legend.position = "none")

ggsave("outputs/Fig2_F_correlation_features_with_lambda.pdf", width=9, height=2.5)
