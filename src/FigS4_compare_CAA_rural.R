suppressPackageStartupMessages({
    library(tidyverse)
    library(limma)
    library(patchwork)
    library(here)
    library(glue)
})
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")


PATHS <- list(
    data = list(
        EXV = "data/EXV_percentLineage.rds",
        PMA = "data/PMA_percentLineage.rds",
        MPL = "data/MPL_percentLineage.rds",
        MET = "data/MET_80ptile.rds",
        GLY = "data/GLY_igg.rds"
    ),
    meta = "data/PMA_metadata.rds",
    CAA = "data/CAA_rural_results.rds",
    results = "outputs/"
)


# load and prepare data ------------------------------------------------------------

# do arcsin-squareroot transformation on cluster proportions
dat.EXV <- asin(sqrt(readRDS(PATHS$data$EXV)/100))
dat.PMA <- asin(sqrt(readRDS(PATHS$data$PMA)/100))
dat.MPL <- asin(sqrt(readRDS(PATHS$data$MPL)/100))
dat.MET <- readRDS(PATHS$data$MET)
dat.GLY <- readRDS(PATHS$data$GLY)


# tidy up names
rownames(dat.EXV) <- str_remove(str_replace_all(rownames(dat.EXV), "\\+", "pos_"), "_$")
rownames(dat.PMA) <- str_remove(str_replace_all(rownames(dat.PMA), "\\+", "pos_"), "_$")
rownames(dat.MPL) <- str_remove(str_replace_all(rownames(dat.MPL), "\\⁺", "pos_"), "_$")

schisto_meta <- readRDS("data/CAA_rural_results.rds")

# EXV -----------------------------------------------------------

dat.sch <- dat.EXV[, colnames(dat.EXV) %in% schisto_meta$sample]

design <- model.matrix(~ CAA_results + age + sex, data = schisto_meta)
fit <- lmFit(dat.sch, design, maxit = 100)
efit <- eBayes(fit, robust = TRUE)
res.exv <- topTable(efit, coef = 2, number = Inf)

p1 <- res.exv |>
    ggplot(aes(logFC, -log10(adj.P.Val))) +
    geom_point() +
    geom_hline(yintercept=-log10(c(0.1)), lty=2) +
    geom_vline(xintercept=0, lty=2)+
    xlim(c(-0.15, 0.15)) +
    ylim(c(0, -log10(0.001))) +
    labs(title="EXV")


# PMA ----------------------------------------------------------------------


dat.pma.sch <- dat.PMA[, colnames(dat.PMA) %in% schisto_meta$sample]
dat.pma.sch <- dat.pma.sch[str_detect(rownames(dat.pma.sch), "[IFN|TNF|IL].pos"),]
dat.pma.sch <- dat.pma.sch[str_detect(rownames(dat.pma.sch), "^[CD4|CD8|NK|B|Mono]"),]

fit <- lmFit(dat.pma.sch, design, method="robust", maxit = 100)
efit <- eBayes(fit, robust = TRUE)
res.pma <- topTable(efit, coef = 2, number = Inf)

p2 <- res.pma |>
    ggplot(aes(logFC, -log10(adj.P.Val))) +
    geom_point() +
    geom_hline(yintercept=-log10(c(0.1)), lty=2) +
    geom_vline(xintercept=0, lty=2)+
    xlim(c(-0.15, 0.15)) +
    ylim(c(0, -log10(0.001))) +
    labs(title="PMA")


# MPL ---------------------------------------------------------------------


dat.mpl.sch <- dat.MPL[, colnames(dat.MPL) %in% schisto_meta$sample]
dat.mpl.sch <- dat.mpl.sch[str_detect(rownames(dat.mpl.sch), "[IFN|TNF|IL].pos"),]
dat.mpl.sch <- dat.mpl.sch[str_detect(rownames(dat.mpl.sch), "^[CD4|CD8|NK|B|Mono]"),]

fit <- lmFit(dat.mpl.sch, design, method="robust", maxit = 100)
efit <- eBayes(fit, robust = TRUE)
res.mpl <- topTable(efit, coef = 2, number = Inf)

p3 <- res.mpl |>
    ggplot(aes(logFC, -log10(adj.P.Val))) +
    geom_point()+
    geom_hline(yintercept=-log10(c(0.1)), lty=2) +
    geom_vline(xintercept=0, lty=2)+
    xlim(c(-0.15, 0.15)) +
    ylim(c(0, -log10(0.001))) +
    labs(title="MPL")

# MET -----------------------------------------------------------------

dat.metflow.sch <- dat.MET[, colnames(dat.MET) %in% schisto_meta$sample]

fit <- lmFit(dat.metflow.sch, design, method="robust", maxit = 1000)
efit.met <- eBayes(fit, robust = TRUE)
res.met <- topTable(efit.met, coef = 2, number = Inf)


p4 <- res.met |>
    rownames_to_column(var="subset") |>
    filter(str_detect(subset, 'CD4T|Mono|NK|B')) |>
    ggplot(aes(logFC, -log10(adj.P.Val))) +
    geom_point()+
    geom_hline(yintercept=-log10(c(0.1)), lty=2) +
    geom_vline(xintercept=0, lty=2)+
    xlim(c(-0.15, 0.15)) +
    ylim(c(0, -log10(0.001))) +
    labs(title="MET")


# GLY -------------------------------------------------------------------


dat.glyco.sch <- dat.GLY[, colnames(dat.GLY) %in% schisto_meta$sample]

fit <- lmFit(dat.glyco.sch, design, method="robust", maxit = 100)
efit <- eBayes(fit, robust = TRUE)
res.gly <- topTable(efit, coef = 2, number = Inf)

p5 <- res.gly |>
    ggplot(aes(logFC, -log10(adj.P.Val))) +
    geom_point()+
    geom_hline(yintercept=-log10(c(0.1)), lty=2) +
    geom_vline(xintercept=0, lty=2)+
    xlim(c(-0.15, 0.15)) +
    ylim(c(0, -log10(0.001)))+
    labs(title="GLY")



# plots ---------------------------------------------------------------------

p1+p2+p3+p4+p5+
    plot_layout(guides="collect", nrow=1) &
    cowplot::theme_half_open() &
    labs(y="-log10(FDR)")

ggsave("outputs/FigS4_CAA_subgroup_volcano.pdf", dpi=600, width=12, height=3)


# TableS6 ----------------------------------------------------

res.exv.top10 <- res.exv |>
    rownames_to_column(var="feature") |>
    select(feature, logFC, t, P.Value, adj.P.Val) |>
    mutate(dataset="EXV", .before="feature") |>
    slice_min(P.Value, n=10)
res.pma.top10 <- res.pma |>
    rownames_to_column(var="feature") |>
    select(feature, logFC, t, P.Value, adj.P.Val) |>
    mutate(dataset="PMA", .before="feature") |>
    slice_min(P.Value, n=10)
res.mpl.top10 <- res.mpl |>
    rownames_to_column(var="feature") |>
    select(feature, logFC, t, P.Value, adj.P.Val) |>
    mutate(dataset="MPL", .before="feature") |>
    slice_min(P.Value, n=10)
res.met.top10 <- res.met |>
    rownames_to_column(var="feature") |>
    select(feature, logFC, t, P.Value, adj.P.Val) |>
    mutate(dataset="MET", .before="feature") |>
    filter(grepl("CD4T|NK|Mono|B", feature)) |>
    slice_min(P.Value, n=10)
res.gly.top10 <- res.gly |>
    rownames_to_column(var="feature") |>
    select(feature, logFC, t, P.Value, adj.P.Val) |>
    mutate(dataset="GLY", .before="feature") |>
    slice_min(P.Value, n=10)

supp.table <- data.table::rbindlist(list(res.exv.top10, res.pma.top10, res.mpl.top10, res.met.top10, res.gly.top10)) |>
    mutate(across(logFC:adj.P.Val, \(x) signif(x, digits=3)))
openxlsx::write.xlsx(supp.table, "outputs/TableS6.xlsx")
