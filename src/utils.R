## PLOTTING
color_lineages <- ggthemes::tableau_color_pal()(8)
color_lineages <- purrr::set_names(color_lineages, c("gdT", "UncT", "CD8T", "CD4T", "B", "Myeloids", "Basophils", "NKILC"))

# color_groups <- c("#0072B2", "#56B4E9", "#7AC36A")
color_groups <- c("#E69F00", "#56B4E9", "#009E73")
color_groups <- purrr::set_names(color_groups, c("RUR", "URB", "EUR"))

# plotScatter <- function(data, x, y, ...) {
#     data %>%
#         ggplot(aes(.data[[x]], .data[[y]], ...)) +
#         geom_scattermore(alpha = .8) +
#         coord_cartesian(xlim = c(0, 8), ylim = c(0, 8))
# }
plotDensity <- function(data, x, y, nbin = 128, pt = 1.1, pal=1) {
    pal <- pal
    color <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))
    dens <- densCols(data[[x]], data[[y]], colramp = color)
    ggplot(data, aes(.data[[x]], .data[[y]])) +
        geom_scattermore(color = dens, pointsize = pt)
}
saveTest <- function(p, width = 6, height = 4, ...) ggsave(filename = paste0(PATHS$plots, "test_plot.png"), plot = p, bg = "white", width = width, height = height, ...)

savePlot <- function(p, name = "test_plot.png", ...) ggsave(filename = paste0(PATHS$plots, name), plot = p, units = "cm", ...)

savePlot2 <- function(p, name = "test_plot", ...){
    cowplot::ggsave2(filename = paste0(PATHS$plots, name, ".png"), plot = p, units = "cm", dpi=600, bg="white", ...)
    cowplot::ggsave2(filename = paste0(PATHS$plots, name, ".pdf"), plot = p, units = "cm", device=cairo_pdf,...)
}


## DIMENSIONALITY REDUCTIONS

runTSNE <- function(data, save_path, overwrite = FALSE) {
    if (file.exists(save_path) & !overwrite) {
        message("Result already exists. Loading...\n")
        return(qs::qread(save_path))
    } else {
        message("Running tSNE...\n")
        tsne <- fftRtsne(as.matrix(data), fast_tsne_path = fast_tsne_path, perplexity_list = c(30, nrow(data)))

        message("Saving results...\n")
        qs::qsave(tsne, save_path)

        message("Done!\n")
        return(tsne)
    }
}

## CLUSTERING
clusterFastPG <- function(data, save_path = NULL, overwrite = FALSE, ...) {
    if (file.exists(save_path) && !overwrite) {
        message("Result already exists. Loading...\n")
        return(qs::qread(save_path))
    } else {
        message("Running FastPG...\n")
        set.seed(42)
        all_knn <- RcppHNSW::hnsw_knn(as.matrix(data), ...)
        ind <- all_knn$idx
        links <- FastPG::rcpp_parallel_jce(ind)
        links <- FastPG::dedup_links(links)
        clusters <- FastPG::parallel_louvain(links)

        message("Saving results...\n")
        qs::qsave(clusters, save_path)

        message("Done!\n")
        return(clusters$communities + 1)
    }
}

clusterSOM <- function(data, save_path, overwrite = FALSE, ...) {
    if (file.exists(save_path) & !overwrite) {
        message("Result already exists. Loading...\n")
        return(qs::qread(save_path))
    } else {
        message("Running FlowSOM...\n")
        set.seed(666)
        som.out <- FlowSOM::SOM(as.matrix(data), ...)
        clusters <- as.factor(som.out$mapping[, 1])

        message("Saving results...\n")
        qs::qsave(clusters, save_path)

        message("Done!\n")
        return(clusters)
    }
}


themeTSNE <- function(font_size = 10, axis_length = 0.3, arrow_length = 0.1, arrow_angle=30, ...) {
    arrow <- grid::arrow(angle=arrow_angle, length = unit(arrow_length, "inches"), type = "closed")
    cowplot::theme_half_open(font_size = font_size, ...) %+replace%
        tidydr::theme_noaxis(
            axis.line.x.bottom = tidydr::element_line2(id = 1, xlength = axis_length, arrow = arrow),
            axis.line.y.left = tidydr::element_line2(id = 2, ylength = axis_length, arrow = arrow),
            axis.title = element_text(hjust = 0.1, vjust=1)
        )
}


# geoms -------------------------------------------------------------------

StatCentSeg <- ggplot2::ggproto("StatCentSeg", Stat,
                                compute_group = function(data, scales, params,
                                                         cfun=median) {
                                    data$xend <- cfun(data$x)
                                    data$yend <- cfun(data$y)
                                    return(data)
                                },
                                required_aes = c("x", "y")
)
stat_centseg <- function(mapping = NULL, data = NULL, geom = "segment",
                         position = "identity", na.rm = FALSE, show.legend = NA,
                         inherit.aes = TRUE, cfun=median, ...) {
    layer(
        stat = StatCentSeg, data = data, mapping = mapping, geom = geom,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, cfun = cfun, ...)
    )
}

# taken from:
# https://gist.github.com/mcgoodman/58c9d1257fd1625954a4ffa1c3301939
pairwise_permanova <- function(sp_matrix, group_var, dist = "bray", adj = "fdr", perm = 10000) {

    require(vegan)

    ## list contrasts
    group_var <- as.character(group_var)
    groups <- as.data.frame(t(combn(unique(group_var), m = 2)))

    contrasts <- data.frame(
        group1 = groups$V1, group2 = groups$V2,
        R2 = NA, F_value = NA, df1 = NA, df2 = NA, p_value = NA
    )

    for (i in seq(nrow(contrasts))) {
        sp_subset <- group_var == contrasts$group1[i] | group_var == contrasts$group2[i]
        contrast_matrix <- sp_matrix[sp_subset,]

        ## fit contrast using adonis
        fit <- vegan::adonis2(
            contrast_matrix ~ group_var[sp_subset],
            method = dist,
            perm = perm
        )

        contrasts$R2[i] <- round(fit$R2[1], digits = 3)
        contrasts$F_value[i] <- round(fit[["F"]][1], digits = 3)
        contrasts$df1[i] <- fit$Df[1]
        contrasts$df2[i] <- fit$Df[2]
        contrasts$p_value[i] <- fit$`Pr(>F)`[1]
    }

    ## adjust p-values for multiple comparisons
    contrasts$p_value <- round(p.adjust(contrasts$p_value, method = adj), digits = 3)

    return(list(
        contrasts = contrasts,
        "p-value adjustment" = adj,
        permutations = perm
    ))
}



get_pvalues <- function(fit, formula){
    require(emmeans)
    ix.sort <- c(1,3,2)

    df <- data.frame(emmeans(fit, formula, adjust="none")$contrasts) %>%
        mutate(p.value=scales::pvalue(p.value, accuracy=0.0001))
    # message(formula)
    if(str_detect(paste0(as.character(eval(pairwise~group2)), collapse=""), "pairwise")) {
        df <- df %>%
            separate(contrast, c("xmin", "xmax"))
        if(nrow(df)/3 == 2) ix.sort <- c(ix.sort, ix.sort+3)
        return(df[ix.sort,])
    } else {
        df
    }
}

geom_signif2 <- function(annotations, margin_top=.03, step_increase=c(0, 0, .05), textsize=7/.pt,
                         extend_line=-.02, tip_length=0.02,
                         comparisons=list(c("RUR", "URB"), c("URB", "EUR"), c("RUR", "EUR")),
                         ...){
    ggpubr::geom_signif(
        annotations=annotations,
        textsize=textsize,
        size = 0.3,
        margin_top=margin_top,
        step_increase=step_increase,
        comparisons=comparisons,
        extend_line=extend_line,
        tip_length=tip_length,
        ...
    )
}


boxplots <- list(
    geom_boxplot(outlier.shape=NA, alpha=.8),
    geom_quasirandom(),
    scale_fill_brewer(type="qual", palette=2),
    labs(x=NULL),
    theme(legend.position = "none")
)



weightedpoly.emmc <- function(levs){
    weights <- c(0.5135135, 0.2702703, 0.2162162)
    width <- 1
    n <- length(levs)
    y <- 1:n - sum(1:n * weights)
    X <- sqrt(weights) * outer(y, seq_len(n) - 1, "^")
    QR <- qr(X)
    z <- QR$qr
    z <- z * (row(z) == col(z))
    raw <- qr.qy(QR, z)/sqrt(weights)
    contr <- sweep(raw, 2L, apply(raw, 2L, function(x) sqrt(sum(x^2))),
                   "/", check.margin = FALSE)
    scores <- seq(1, width * n, by = width)
    scores <- scores - sum(scores * weights)
    contr[, 2] <- contr[, 2] * sqrt(sum(scores^2))
    dn <- paste0("^", 1L:n - 1L)
    dn[2:min(4, n)] <- c(".L", ".Q", ".C")[1:min(3, n - 1)]
    colnames(contr) <- dn
    data.frame(contr[, -1, drop = FALSE])
}
