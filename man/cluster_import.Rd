% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_import.R
\name{cluster_import}
\alias{cluster_import}
\title{Import cluster results}
\usage{
cluster_import(
  spe,
  cluster_dir = file.path(tempdir(), "exported_clusters"),
  prefix = "imported_",
  overwrite = TRUE
)
}
\arguments{
\item{spe}{A
\link[SpatialExperiment:SpatialExperiment]{SpatialExperiment-class}
object. See \code{\link[=fetch_data]{fetch_data()}} for how to download some example objects or
\code{\link[=read10xVisiumWrapper]{read10xVisiumWrapper()}} to read in \code{spaceranger --count} output files and
build your own \code{spe} object.}

\item{cluster_dir}{A \code{character(1)} specifying the output directory, similar
to the \code{outs/analysis/clustering} produced by SpaceRanger.}

\item{prefix}{A \code{character(1)} specifying the prefix to use when naming
these new cluster variables.}

\item{overwrite}{A \code{logical(1)} indicating whether to overwrite the \code{spe$key}.}
}
\value{
A
\link[SpatialExperiment:SpatialExperiment]{SpatialExperiment-class} object
with the imported clusters appended on the \code{colData()}.
}
\description{
This function imports previously exported clustering results with
\code{cluster_export()} and adds them to the \code{colData()} slot of your
\link[SpatialExperiment:SpatialExperiment]{SpatialExperiment-class}
object.
}
\examples{
if (enough_ram()) {
    ## Obtain the necessary data
    if (!exists("spe")) spe <- fetch_data("spe")

    ## Export two cluster variables
    cluster_export(spe, "spatialLIBD")
    cluster_export(spe, "GraphBased")

    ## Re-import them
    colData(cluster_import(spe))
}
}
\seealso{
Other cluster export/import utility functions: 
\code{\link{cluster_export}()}
}
\concept{cluster export/import utility functions}
