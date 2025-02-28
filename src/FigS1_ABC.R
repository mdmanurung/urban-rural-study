
library(tidyverse)

PATHS <- list(
    diffcounts = "data/diffcounts.rds",
    GLY = "data/GLY_igg.rds",
    CAA = "data/CAA_rural_results.rds",
    luminex = "data/data_luminex.rds"
)
noLegend <- list(theme(legend.position = "none"))


# FigS1A ------------------------------------------------------------------

diff <- readRDS(PATHS$diffcounts)
p <- diff %>%
    ggplot(aes(group3, `Eosino_%`)) +
    geom_boxplot(outlier.shape=NA, aes(fill=group3), alpha=.5) +
    geom_quasirandom(aes(shape=group3)) +
    # compare_groups_pairwise() +
    labs(x=NULL, y="% Eosinophils") +
    scale_fill_groups() +
    coord_cartesian(ylim=c(0,33)) +
    noLegend
ggsave("outputs/FigS1_A.pdf")



# FigS1B ------------------------------------------------------------------

caa <- readRDS(PATHS$CAA)
GLY <- readRDS(PATHS$GLY)

df.GLY <- GLY |>
    t() |>
    as.data.frame() |>
    rownames_to_column(var="sample") |>
    left_join(caa) |>
    drop_na(CAA_results)

# adjust for age and sex
fit <- lm(sin(IgGI_Galactosylation)^2*100~CAA_results+age+sex, data=df.GLY)
emm <- emmeans(fit, pairwise ~ CAA_results)

p <- df.GLY |>
    ggplot(aes(CAA_results, sin(IgGI_Galactosylation)^2*100)) +
    ggprism::add_pvalue(data.frame(group1="CAA+", group2="CAA-", p="0.21", y.position=68), tip.length=0.01) +
    geom_errorbar(data = data.frame(emm$emmeans), aes(y=emmean, ymin=lower.CL, ymax=upper.CL), width=.5) +
    geom_point(data = data.frame(emm$emmeans), aes(y=emmean), size=4, color="black") +
    ggbeeswarm::geom_quasirandom() +
    labs(x="", y="% IgG1 Fc-galactosylation") +
    scale_x_discrete(labels=c("CAA (+)", "CAA (-)")) +
    cowplot::theme_cowplot(font_size=12)
ggsave("outputs/FigS1_B.pdf", p, dpi=600, width=2, height=3)



# FigS1_C -----------------------------------------------------------------

# replace with path to 'data_luminex.tsv'
data <- readRDS(PATHS$luminex)


# calculate trend test pvalues
pvalues_jt <- sapply(colnames(data)[4:12], function(x){
    message(x)
    set.seed(42)
    pvals <- clinfun::jonckheere.test(data[[x]], ordered(data$group), nperm=2000)
    pvals$p.value
})
fdr_jt <- p.adjust(pvalues_jt)
(sig_markers <- names(sort(fdr_jt[fdr_jt < .1])))


# calculate pairwise p values
pvalues_kw <- lapply(colnames(data)[4:12], function(x){
    message(x)
    set.seed(42)
    pvals <- rstatix::pairwise_wilcox_test(data=data, as.formula(paste0(x, " ~ group")), p.adjust.method="none")
    pvals |>
        mutate(comparisons = paste0("wilcoxPval_", group1, "vs.", group2)) |>
        select(comparisons, p.adj) |>
        mutate(comparisons = str_replace_all(comparisons, "\n", "_")) |>
        pivot_wider(names_from=comparisons, values_from=p.adj)
})
pvalues_kw <- do.call(rbind, pvalues_kw)


group_means <- data |>
    group_by(group) |>
    summarise(across(ccl24_eotaxin2_mpif2:cd40ligand_tnfsf5, mean)) |>
    mutate(group = str_replace_all(group, "\n", "_"),
           group = paste0("mean_", group)) |>
    column_to_rownames(var="group") |>
    t() |>
    as.data.frame() |>
    rownames_to_column(var="analyte")
group_means$trend_test_FDR <- fdr_jt
group_means <- cbind(group_means, pvalues_kw)
group_means <- group_means |>
    arrange(trend_test_FDR)

# supplementary table 4
write.xlsx(group_means, "outputs/TableS3_luminex_stats.xlsx")


sort(fdr_jt[fdr_jt < .1])
# complement_component_c5a  il8_cxcl8
# 0.045                     0.064
p <- data |>
    pivot_longer(ccl24_eotaxin2_mpif2:cd40ligand_tnfsf5) |>
    filter(name %in% c("complement_component_c5a", "il8_cxcl8", "ccl2_je_mcp_1")) |>
    mutate(name = factor(name, levels=c("complement_component_c5a", "il8_cxcl8", "ccl2_je_mcp_1"))) |>
    ggplot(aes(group, value)) +
    geom_boxplot(outlier.shape=NA, aes(fill=group), alpha=.4) +
    ggbeeswarm::geom_quasirandom(aes(color=group)) +
    facet_wrap(~name, scales="free_y", nrow=1) +
    ggpubr::stat_pwc(hide.ns=F, method="wilcox_test", label.size=8/.pt, step.increase=.08, tip.length=.01) +
    cowplot::theme_half_open(10) +
    theme(strip.background = element_blank(), legend.position = "none") +
    labs(x=NULL, y="Log2 Conc. (pg/mL)") +
    scale_color_manual(values=ggthemes::colorblind_pal()(8)[c(2,3,4)]) +
    scale_fill_manual(values=ggthemes::colorblind_pal()(8)[c(2,3,4)])
ggsave("outputs/FigS1_C.pdf", width=6, height=2.5)

