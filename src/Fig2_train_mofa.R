suppressPackageStartupMessages({
    library(MOFA2)
    library(tidyverse)
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
    results = "outputs/",
    models = "outputs/mofa/trained_models/"
)

# for scaling
tst <- function(matrix) t(scale(t(matrix)))

## PARAMETERS
prop_top_features <- 1
mofa_iter <- 50

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

meta <- readRDS(PATHS$meta) %>%
    rename(sample = patient_id) %>%
    as.data.frame() %>%
    mutate(group3 = factor(group3, levels=c("RUR", "URB", "EUR")))
rownames(meta) <- meta$sample


# prepare data to correct format for MOFA
# ensure same order of samples across tables


# exclude clusters not enriched for cytokines for PMA and MPL data
ix <- str_which(rownames(dat.PMA), "IFN|IL4|TNF|IL1")
dat.PMA <- dat.PMA[ix, colnames(dat.EXV)]

ix <- str_which(rownames(dat.MPL), "IFN|IL4|IL1|TNF")
dat.MPL <- dat.MPL[ix, colnames(dat.EXV)]

ix <- !grepl("GLUT1", rownames(dat.MET))
dat.MET <- dat.MET[ix, colnames(dat.EXV)]

dat.GLY <- dat.GLY[, colnames(dat.EXV)]


datList <- list(
    EXV = as.matrix(tst(dat.EXV)),
    PMA = as.matrix(tst(dat.PMA)),
    MPL = as.matrix(tst(dat.MPL)),
    MET = as.matrix(tst(dat.MET)),
    GLY = as.matrix(tst(dat.GLY))
)



# train MOFA model --------------------------------------------------------


(MOFAobj <- create_mofa(datList))

(data_opts <- get_default_data_options(MOFAobj))
data_opts$scale_views <- TRUE

(model_opts <- get_default_model_options(MOFAobj))
model_opts$num_factors <- 9
model_opts$ard_factors <- TRUE
model_opts$spikeslab_weights <- TRUE

path_mofa_models <- "outputs/mofa/mofa_models.rds"
if (file.exists(path_mofa_models)) {

    mofa_models <- readRDS(path_mofa_models)

} else {

    mofa_models <- lapply(1:mofa_iter, function(i){
        (train_opts <- get_default_training_options(MOFAobj))
        train_opts$convergence_mode <- "slow"
        train_opts$seed <- runif(1, 0, 99999)
        train_opts$maxiter <- 2000
        train_opts$drop_factor_threshold <- 0.01

        MOFAobj <- prepare_mofa(
            object = MOFAobj,
            data_options = data_opts,
            model_options = model_opts,
            training_options = train_opts
        )
        samples_metadata(MOFAobj) <- meta

        outfile = here(PATHS$models, glue("mofa_trained_{i}.hdf5"))
        MOFAobj.trained <- run_mofa(MOFAobj, outfile, use_basilisk=TRUE)
    })

}


# extract best model and latent factors -----------------------------------

p <- MOFA2::compare_elbo(mofa_models)
p$data |>
    mutate(ELBO = ELBO/max(ELBO)) |>
    mutate(model = fct_reorder(model, ELBO, .desc=TRUE)) |>
    ggplot(aes(model, ELBO)) +
    geom_col() +
    coord_cartesian(ylim=c(.98, 1)) +
    ggpubr::rotate_x_text()

# pick model with the highest elbo, which is model 8
selected_model <- MOFA2::select_model(models)
saveRDS(selected_model, "outputs/mofa/final_mofa_model.rds")


# extract MOFA latent factors for downstream analyses
latent.factors.coords <- get_factors(mofa_models[[8]], scale=FALSE)$group1
latent.factors <- latent.factors.coords |>
    as.data.frame() |>
    rownames_to_column(var="sample") |>
    left_join(meta, by="sample")
saveRDS(latent.factors, "outputs/mofa/final_mofa_model_factors.rds")



# plot MOFA schematics ----------------------------------------------------


mofa <- readRDS("outputs/mofa/final_mofa_model.rds")
plots <- MOFA2::plot_variance_explained(mofa, plot_total=TRUE, max_r2=10)
p1 <- plots[[1]]
p2 <- plots[[2]]

plot1 <- p1$data |>
    mutate(view = fct_rev(view)) |>
    mutate(factor = str_replace(factor, "Factor", "F")) |>

    ggplot(aes(factor, y=view)) +
    geom_tile(aes(fill=value), color="gray90") +
    scale_fill_viridis_c(name="Var. Explained (%)",option="G", limits=c(0,10)) +
    coord_cartesian(expand=FALSE) +
    labs(x="Latent Factors", y=NULL, subtitle=NULL)+
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.position = "top") +
    guides(colour=guide_colourbar(direction="horizontal", title.position="top"))

plot2 <- p2$data |>
    mutate(view = fct_rev(view)) |>

    ggplot(aes(R2, view)) +
    geom_col(fill="gray70") +
    coord_cartesian(expand=FALSE) +
    labs(x="Total Var. Explained", y=NULL, subtitle=NULL)+
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

plot1+plot2+plot_layout(widths=c(2,1))
ggsave("outputs/Fig2_B_mofaVariance.pdf", dpi=600, width=3, height=3)

