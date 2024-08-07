% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_app.R
\name{run_app}
\alias{run_app}
\title{Run the spatialLIBD Shiny Application}
\usage{
run_app(
  spe = fetch_data(type = "spe"),
  sce_layer = fetch_data(type = "sce_layer"),
  modeling_results = fetch_data(type = "modeling_results"),
  sig_genes = sig_genes_extract_all(n = nrow(sce_layer), modeling_results =
    modeling_results, sce_layer = sce_layer),
  docs_path = system.file("app", "www", package = "spatialLIBD"),
  title = "spatialLIBD",
  spe_discrete_vars = c("spatialLIBD", "GraphBased", "ManualAnnotation", "Maynard",
    "Martinowich", paste0("SNN_k50_k", 4:28), "SpatialDE_PCA", "SpatialDE_pool_PCA",
    "HVG_PCA", "pseudobulk_PCA", "markers_PCA", "SpatialDE_UMAP", "SpatialDE_pool_UMAP",
    "HVG_UMAP", "pseudobulk_UMAP", "markers_UMAP", "SpatialDE_PCA_spatial",
    "SpatialDE_pool_PCA_spatial", "HVG_PCA_spatial", "pseudobulk_PCA_spatial",
    "markers_PCA_spatial", "SpatialDE_UMAP_spatial", "SpatialDE_pool_UMAP_spatial",
    "HVG_UMAP_spatial", "pseudobulk_UMAP_spatial", 
     "markers_UMAP_spatial"),
  spe_continuous_vars = c("cell_count", "sum_umi", "sum_gene", "expr_chrM",
    "expr_chrM_ratio"),
  default_cluster = "spatialLIBD",
  auto_crop_default = TRUE,
  is_stitched = FALSE,
  ...
)
}
\arguments{
\item{spe}{Defaults to the output of
\code{fetch_data(type = 'spe')}. This is a
\link[SpatialExperiment:SpatialExperiment]{SpatialExperiment-class}
object with the spot-level Visium data and information required for
visualizing the histology. See \code{\link[=fetch_data]{fetch_data()}} for more details.}

\item{sce_layer}{Defaults to the output of
\code{fetch_data(type = 'sce_layer')}. This is a
\linkS4class{SingleCellExperiment}
object with the spot-level Visium data compressed via pseudo-bulking to the
layer-level (group-level) resolution. See \code{\link[=fetch_data]{fetch_data()}} for more details.}

\item{modeling_results}{Defaults to the output of
\code{fetch_data(type = 'modeling_results')}. This is a list of tables with the
columns \verb{f_stat_*} or \verb{t_stat_*} as well as \verb{p_value_*} and \verb{fdr_*} plus
\code{ensembl}. The column name is used to extract the statistic results, the
p-values, and the FDR adjusted p-values. Then the \code{ensembl} column is used
for matching in some cases. See \code{\link[=fetch_data]{fetch_data()}} for more details.}

\item{sig_genes}{The output of \code{\link[=sig_genes_extract_all]{sig_genes_extract_all()}} which is a table
in long format with the modeling results. You can subset this if the object
requires too much memory.}

\item{docs_path}{A \code{character(1)} specifying the path to the directory
containing the website documentation files. The directory has to contain
the files: \code{documentation_sce_layer.md}, \code{documentation_spe.md},
\code{favicon.ico}, \code{footer.html} and \code{README.md}.}

\item{title}{A character(1) specifying the title for the app.}

\item{spe_discrete_vars}{A \code{character()} vector of discrete variables that
will be available to visualize in the app. Basically, the set of variables
with spot-level groups. They will have to be present in \code{colData(spe)}.}

\item{spe_continuous_vars}{A \code{character()} vector of continuous variables
that will be available to visualize in the app using the same scale
as genes. They will have to be present in \code{colData(sce)}.}

\item{default_cluster}{A \code{character(1)} with the name of the main cluster
(discrete) variable to use. It will have to be present in both \code{colData(spe)}
and \code{colData(sce_layer)}.}

\item{auto_crop_default}{A \code{logical(1)} specifying the default value for
automatically cropping the images. Set this to \code{FALSE} if your images do not
follow the Visium grid size expectations, which are key for enabling
auto-cropping.}

\item{is_stitched}{A \code{logical(1)} vector: If \code{TRUE}, expects a
\link[SpatialExperiment:SpatialExperiment]{SpatialExperiment-class} built
with \code{visiumStitched::build_spe()}.
\url{http://research.libd.org/visiumStitched/reference/build_spe.html}; in
particular, expects a logical colData column \code{exclude_overlapping}
specifying which spots to exclude from the plot. Sets \code{auto_crop = FALSE}.}

\item{...}{Other arguments passed to the list of golem options for running
the application.}
}
\value{
A \link[shiny:shiny.appobj]{shiny.appobj} that contains the input data.
}
\description{
This function runs the shiny application that allows users to interact
with the Visium spatial transcriptomics data from LIBD (by default) or
any other data that you have shaped according to our object structure.
}
\details{
If you don't have the pseudo-bulked analysis results like we computed them
in our project \url{https://doi.org/10.1038/s41593-020-00787-0} you can
set \code{sce_layer}, \code{modeling_results} and \code{sig_genes} to \code{NULL}. Doing so
will disable the pseudo-bulked portion of the web application. See the
examples for one such case as well as the vignette that describes how
you can use \code{spatialLIBD} with public data sets provided by 10x Genomics.
That vignette is available at
\url{http://research.libd.org/spatialLIBD/articles/TenX_data_download.html}.
}
\examples{
\dontrun{
## The default arguments will download the data from the web
## using fetch_data(). If this is the first time you have run this,
## the files will need to be cached by ExperimentHub. Otherwise it
## will re-use the files you have previously downloaded.
if (enough_ram(4e9)) {
    ## Obtain the necessary data
    if (!exists("spe")) spe <- fetch_data("spe")

    ## Create the interactive website
    run_app(spe)

    ## You can also run a custom version without the pseudo-bulked
    ## layer information. This is useful if you are only interested
    ## in the spatial transcriptomics features.
    run_app(spe,
        sce_layer = NULL, modeling_results = NULL, sig_genes = NULL,
        title = "spatialLIBD without layer info"
    )

    ## When using shinyapps.io aim for less than 3 GB of RAM with your
    ## objects. Check each input object with:
    ## lobstr::obj_size(x)
    ## Do not create the large input objects on the app.R script before
    ## subsetting them. Do this outside app.R since the app.R script is
    ## run at shinyapps.io, so subsetting on that script to reduce the
    ## memory load is pointless. You have to do it outside of app.R.
}

## How to run locally the spatialDLPFC Sp09 spatialLIBD app. That is,
## from http://research.libd.org/spatialDLPFC/#interactive-websites
## how to run https://libd.shinyapps.io/spatialDLPFC_Visium_Sp09 locally.
if (enough_ram(9e9)) {
    ## Download the 3 main objects needed
    spe <- fetch_data("spatialDLPFC_Visium")
    sce_pseudo <- fetch_data("spatialDLPFC_Visium_pseudobulk")
    modeling_results <- fetch_data("spatialDLPFC_Visium_modeling_results")

    ## These are optional commands to further reduce the memory required.
    #
    ## Keep only the "lowres" images. Reduces the object from 6.97 GB to 4.59 GB
    # imgData(spe) <- imgData(spe)[imgData(spe)$image_id == "lowres", ]
    ## Drop the regular counts (keep only the logcounts). Reduces the object
    ## from 4.59 GB to 2.45 GB.
    # counts(spe) <- NULL

    ## For sig_genes_extract_all() to work
    sce_pseudo$spatialLIBD <- sce_pseudo$BayesSpace
    ## Compute the significant genes
    sig_genes <- sig_genes_extract_all(
        n = nrow(sce_pseudo),
        modeling_results = modeling_results,
        sce_layer = sce_pseudo
    )
    ## Reduce the memory from 423.73 MB to 78.88 MB
    lobstr::obj_size(sig_genes)
    sig_genes$in_rows <- NULL
    sig_genes$in_rows_top20 <- NULL
    lobstr::obj_size(sig_genes)

    ## Specify the default variable
    spe$BayesSpace <- spe$BayesSpace_harmony_09
    ## Get all variables
    vars <- colnames(colData(spe))

    ## Set default cluster colors
    colors_BayesSpace <- Polychrome::palette36.colors(28)
    names(colors_BayesSpace) <- c(1:28)
    m <- match(as.character(spe$BayesSpace_harmony_09), names(colors_BayesSpace))
    stopifnot(all(!is.na(m)))
    spe$BayesSpace_colors <- spe$BayesSpace_harmony_09_colors <- colors_BayesSpace[m]

    ## Download documentation files we use
    temp_www <- file.path(tempdir(), "www")
    dir.create(temp_www)
    download.file(
        "https://raw.githubusercontent.com/LieberInstitute/spatialDLPFC/main/README.md",
        file.path(temp_www, "README.md")
    )
    download.file(
        "https://raw.githubusercontent.com/LieberInstitute/spatialDLPFC/main/code/deploy_app_k09/www/documentation_sce_layer.md",
        file.path(temp_www, "documentation_sce_layer.md")
    )
    download.file(
        "https://raw.githubusercontent.com/LieberInstitute/spatialDLPFC/main/code/deploy_app_k09/www/documentation_spe.md",
        file.path(temp_www, "documentation_spe.md")
    )
    download.file(
        "https://raw.githubusercontent.com/LieberInstitute/spatialDLPFC/main/img/favicon.ico",
        file.path(temp_www, "favicon.ico")
    )
    download.file(
        "https://raw.githubusercontent.com/LieberInstitute/spatialDLPFC/main/code/deploy_app_k09/www/footer.html",
        file.path(temp_www, "footer.html")
    )
    list.files(temp_www)

    ## Run the app locally
    run_app(
        spe,
        sce_layer = sce_pseudo,
        modeling_results = modeling_results,
        sig_genes = sig_genes,
        title = "spatialDLPFC, Visium, Sp09",
        spe_discrete_vars = c( # this is the variables for the spe object not the sce_pseudo object
            "BayesSpace",
            "ManualAnnotation",
            vars[grep("^SpaceRanger_|^scran_", vars)],
            vars[grep("^BayesSpace_harmony", vars)],
            vars[grep("^BayesSpace_pca", vars)],
            "graph_based_PCA_within",
            "PCA_SNN_k10_k7",
            "Harmony_SNN_k10_k7",
            "manual_layer_label",
            "wrinkle_type",
            "BayesSpace_colors"
        ),
        spe_continuous_vars = c(
            "sum_umi",
            "sum_gene",
            "expr_chrM",
            "expr_chrM_ratio",
            vars[grep("^VistoSeg_", vars)],
            vars[grep("^layer_", vars)],
            vars[grep("^broad_", vars)]
        ),
        default_cluster = "BayesSpace",
        docs_path = temp_www
    )
}
## See also:
## * https://github.com/LieberInstitute/spatialDLPFC/tree/main/code/deploy_app_k09
## * https://github.com/LieberInstitute/spatialDLPFC/tree/main/code/deploy_app_k09_position
## * https://github.com/LieberInstitute/spatialDLPFC/tree/main/code/deploy_app_k09_position_noWM
## * https://github.com/LieberInstitute/spatialDLPFC/tree/main/code/deploy_app_k16
## * https://github.com/LieberInstitute/spatialDLPFC/tree/main/code/analysis_IF/03_spatialLIBD_app


## Example for an object with multiple capture areas stitched together with
## <http://research.libd.org/visiumStitched/>.
spe_stitched <- fetch_data("Visium_LS_spe")

## Inspect this object
spe_stitched

## Notice the use of "exclude_overlapping"
table(spe_stitched$exclude_overlapping, useNA = "ifany")

## Run the app with this stitched data
run_app(
    spe = spe_stitched,
    sce_layer = NULL, modeling_results = NULL, sig_genes = NULL,
    title = "visiumStitched example data",
    spe_discrete_vars = c("capture_area", "scran_quick_cluster", "ManualAnnotation"),
    spe_continuous_vars = c("sum_umi", "sum_gene", "expr_chrM", "expr_chrM_ratio"),
    default_cluster = "scran_quick_cluster",
    is_stitched = TRUE
)
}
}
