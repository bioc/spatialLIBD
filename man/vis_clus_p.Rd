% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis_clus_p.R
\name{vis_clus_p}
\alias{vis_clus_p}
\title{Sample spatial cluster visualization workhorse function}
\usage{
vis_clus_p(
  spe,
  d,
  clustervar,
  sampleid = unique(spe$sample_id)[1],
  colors,
  spatial,
  title,
  image_id = "lowres",
  alpha = NA,
  point_size = 2,
  auto_crop = TRUE,
  na_color = "#CCCCCC40"
)
}
\arguments{
\item{spe}{A
\link[SpatialExperiment:SpatialExperiment]{SpatialExperiment-class}
object. See \code{\link[=fetch_data]{fetch_data()}} for how to download some example objects or
\code{\link[=read10xVisiumWrapper]{read10xVisiumWrapper()}} to read in \code{spaceranger --count} output files and
build your own \code{spe} object.}

\item{d}{A \code{data.frame()} with the sample-level information. This is
typically obtained using \code{cbind(colData(spe), spatialCoords(spe))}.}

\item{clustervar}{A \code{character(1)} with the name of the \code{colData(spe)}
column that has the cluster values.}

\item{sampleid}{A \code{character(1)} specifying which sample to plot from
\code{colData(spe)$sample_id} (formerly \code{colData(spe)$sample_name}).}

\item{colors}{A vector of colors to use for visualizing the clusters
from \code{clustervar}. If the vector has names, then those should match the
values of \code{clustervar}.}

\item{spatial}{A \code{logical(1)} indicating whether to include the histology
layer from \code{\link[=geom_spatial]{geom_spatial()}}. If you plan to use
\link[plotly:ggplotly]{ggplotly()} then it's best to set this to \code{FALSE}.}

\item{title}{The title for the plot.}

\item{image_id}{A \code{character(1)} with the name of the image ID you want to
use in the background.}

\item{alpha}{A \code{numeric(1)} in the \verb{[0, 1]} range that specifies the
transparency level of the data on the spots.}

\item{point_size}{A \code{numeric(1)} specifying the size of the points. Defaults
to \code{1.25}. Some colors look better if you use \code{2} for instance.}

\item{auto_crop}{A \code{logical(1)} indicating whether to automatically crop
the image / plotting area, which is useful if the Visium capture area is
not centered on the image and if the image is not a square.}

\item{na_color}{A \code{character(1)} specifying a color for the NA values.
If you set \code{alpha = NA} then it's best to set \code{na_color} to a color that has
alpha blending already, which will make non-NA values pop up more and the NA
values will show with a lighter color. This behavior is lost when \code{alpha} is
set to a non-\code{NA} value.}
}
\value{
A \link[ggplot2:ggplot]{ggplot2} object.
}
\description{
This function visualizes the clusters for one given sample at the spot-level
using (by default) the histology information on the background. This is the
function that does all the plotting behind \code{\link[=vis_clus]{vis_clus()}}. To visualize
gene-level (or any continuous variable) use \code{\link[=vis_gene_p]{vis_gene_p()}}.
}
\examples{

if (enough_ram()) {
    ## Obtain the necessary data
    if (!exists("spe")) spe <- fetch_data("spe")
    spe_sub <- spe[, spe$sample_id == "151673"]

    ## Use the manual color palette by Lukas M Weber
    ## Don't plot the histology information
    p <- vis_clus_p(
        spe = spe_sub,
        d = as.data.frame(cbind(colData(spe_sub), SpatialExperiment::spatialCoords(spe_sub)), optional = TRUE),
        clustervar = "layer_guess_reordered",
        sampleid = "151673",
        colors = libd_layer_colors,
        title = "151673 LIBD Layers",
        spatial = FALSE
    )
    print(p)

    ## Clean up
    rm(spe_sub)
}
}
\seealso{
Other Spatial cluster visualization functions: 
\code{\link{frame_limits}()},
\code{\link{vis_clus}()},
\code{\link{vis_grid_clus}()}
}
\concept{Spatial cluster visualization functions}
