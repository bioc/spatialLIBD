% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read10xVisiumAnalysis.R
\name{read10xVisiumAnalysis}
\alias{read10xVisiumAnalysis}
\title{Load analysis data from a 10x Genomics Visium experiment}
\usage{
read10xVisiumAnalysis(
  samples = "",
  sample_id = paste0("sample", sprintf("\%02d", seq_along(samples)))
)
}
\arguments{
\item{samples}{Passed to \code{\link[SpatialExperiment:read10xVisium]{SpatialExperiment::read10xVisium()}}.}

\item{sample_id}{Passed to \code{\link[SpatialExperiment:read10xVisium]{SpatialExperiment::read10xVisium()}}.}
}
\value{
A named \code{list()} with the information about the clustering and the
dimension reduction (projections) from the SpaceRanger output by 10x
Genomics.
}
\description{
This function expands \code{\link[SpatialExperiment:read10xVisium]{SpatialExperiment::read10xVisium()}} by reading
analysis outputs from SpaceRanger by 10x Genomics.
}
\details{
You might want to use \code{read10xVisiumWrapper()} instead of using this
function directly.
}
\examples{
## See 'Using spatialLIBD with 10x Genomics public datasets' for
## a full example using this function.
if (interactive()) {
    browseVignettes(package = "spatialLIBD")
}

## Note that ?SpatialExperiment::read10xVisium doesn't include all the files
## we need to illustrate read10xVisiumWrapper().
}
\seealso{
Other Utility functions for reading data from SpaceRanger output by 10x
Genomics: 
\code{\link{add10xVisiumAnalysis}()},
\code{\link{read10xVisiumWrapper}()}
}
\concept{Utility functions for reading data from SpaceRanger output by 10x
Genomics}
