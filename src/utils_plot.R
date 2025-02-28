
# shortcut for pairwise test of group means
compare_groups_pairwise <- function(test="t.test", size=8/.pt, ...){
    ggpubr::stat_compare_means(
        comparisons = list(c("RUR", "URB"), c("URB", "EUR"), c("RUR", "EUR")),
        method = test,
        size = size,
        ...
    )
}


scale_color_groups <- function(...) scale_color_manual(values = color_groups, ...)
scale_fill_groups <- function(...) scale_fill_manual(values = color_groups, ...)

scale_color_expression <- function(limits=c(0,5), ...) scale_color_distiller(palette="RdYlBu", limits=limits, oob=scales::oob_squish, ...)


#' Plot boxplots with pairwise tests and trend test
#'
#' @param data
#' @param y
#' @param test
#' @param trend_test
#'
#' @return a ggplot2 boxplot
#' @export
#'
#' @examples
plotBoxplots <- function(dat, y, test="t.test"){
    stopifnot("group3" %in% colnames(dat))
    stopifnot(y %in% colnames(dat))
    p <- dat %>%
        ggplot(aes(group3, .data[[y]], fill=group3)) +
        geom_boxplot(outlier.shape=NA, alpha=.6) +
        ggbeeswarm::geom_quasirandom(aes(shape=group3)) +
        scale_fill_groups() +
        labs(x=NULL)

    if (!is.null(test)) {
        p <- p + compare_groups_pairwise(test=test)
    }

    p
}

plotScatter <- function(dat, x, y, centroid=F, regress=F) {
    stopifnot(c("group3", x, y) %in% colnames(dat))
    p <- dat %>%
        ggplot(aes(.data[[x]], .data[[y]], color=group3)) +
        geom_point() +
        scale_color_groups()

    if (centroid) p <- p + stat_centseg()

    # if(regress) p <- p + geom_smooth(method="lm", se=F)
    if(regress) p <- p + geom_smooth(method="lmrob", se=F)

    p
}

save_heatmap_pdf <- function(x, filename, width=7, height=7, ...) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    cairo_pdf(filename, width=width, height=height, ...)
    ComplexHeatmap::draw(x)
    dev.off()
}

tst <- function(x) t(scale(t(x)))
