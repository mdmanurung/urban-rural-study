

library(tidyverse)
library(here)
library(openxlsx)
library(ggpubr)
library(patchwork)
library(ggbeeswarm)

source("src/utils.R")
source("src/utils_plot.R")
theme_set(theme_pubr(8, border = TRUE))

PATHS <- list(
    meta = "data/PMA_metadata.rds",
    GLY = "data/GLY_igg.rds",
    plots = "outputs/"
)

meta <- readRDS(PATHS$meta) %>%
    mutate(group3 = factor(group3, levels=c("RUR", "URB", "EUR")))


# Fig1B -------------------------------------------------------------------

# age distribution
p1 <- meta %>%
    ggplot(aes(group3, age)) +
    geom_boxplot(outlier.shape=NA, aes(fill=group3), alpha=.5) +
    geom_quasirandom(aes(shape=group3)) +
    labs(x=NULL, y="Age (years)", subtitle = "Kruskal-Wallis P = 0.98") +
    stat_compare_means(size=6/.pt) +
    theme(legend.position = "none") +
    scale_fill_groups()


chisq.test(with(meta, table(group3, sex)))
# X-squared = 0.31144, df = 2, p-value = 0.8558


p2 <- meta %>%
    count(group3, sex) %>%
    ggplot(aes(group3, n, fill=sex)) +
    geom_bar(position="fill", stat="identity") +
    geom_text(aes(label = n), position=position_fill(vjust=.5), color="white", size=4)+
    theme(legend.position = "right") +
    labs(x = NULL, y="Proportion", subtitle="Chisq P = 0.86") +
    coord_cartesian(expand=F)

p <- p1+p2
savePlot(p, "Fig1_B_AgeSex.pdf", width=10, height=5)



# Fig1C -------------------------------------------------------------------

GLY <- readRDS(PATHS$GLY)
df.GLY <- GLY |>
    t() |>
    as.data.frame() |>
    rownames_to_column(var="patient_id") |>
    left_join(meta)

# adjust for age and sex
fit <- gls(IgGI_Galactosylation~group3+age+sex,
           weights = varIdent(form = ~ group3),
           data=df.GLY |> drop_na(IgGI_Galactosylation))
emm.trends <- emmeans(fit, poly ~ group3)
emm.trends$contrasts

emm <- emmeans(fit, pairwise ~ group3)
# emm$contrasts
# contrast  estimate     SE df t.ratio p.value
# RUR - URB  -0.0378 0.0264 24  -1.431  0.3415
# RUR - EUR  -0.0977 0.0323 24  -3.020  0.0157 <---
# URB - EUR  -0.0599 0.0378 24  -1.584  0.2717


p <- df.GLY |>
    ggplot(aes(group3, IgGI_Galactosylation)) +
    geom_boxplot(outlier.shape=NA, aes(fill=group3), alpha=.8) +
    geom_quasirandom() +
    ggprism::add_pvalue(data.frame(group1="RUR", group2="EUR", p="0.0157", y.position=.98), tip.length=0.01) +
    labs(x=NULL, y="% IgG1 Galactosylation") +
    scale_fill_groups()
savePlot(p, "Fig1_C_IgG1Glyco.pdf", width=5, height=5)
