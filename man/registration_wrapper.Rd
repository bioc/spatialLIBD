% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/registration_wrapper.R
\name{registration_wrapper}
\alias{registration_wrapper}
\title{Spatial registration: wrapper function}
\usage{
registration_wrapper(
  sce,
  var_registration,
  var_sample_id,
  covars = NULL,
  gene_ensembl = NULL,
  gene_name = NULL,
  prefix = ""
)
}
\arguments{
\item{sce}{A
\link[SingleCellExperiment:SingleCellExperiment]{SingleCellExperiment-class}
object or one that inherits its properties.}

\item{var_registration}{A \code{character(1)} specifying the \code{colData(sce)}
variable of interest against which will be used for computing the relevant
statistics.}

\item{var_sample_id}{A \code{character(1)} specifying the \code{colData(sce)} variable
with the sample ID.}

\item{covars}{A \code{character()} with names of sample-level covariates.}

\item{gene_ensembl}{A \code{character(1)} specifying the \code{rowData(sce_pseudo)}
column with the ENSEMBL gene IDs. This will be used by \code{layer_stat_cor()}.}

\item{gene_name}{A \code{character(1)} specifying the \code{rowData(sce_pseudo)}
column with the gene names (symbols).}

\item{prefix}{A \code{character(1)} specifying the prefix to use for the
F-statistics column. This is particularly useful if you will run this
function more than once and want to be able to merge the results.}
}
\value{
A \code{list()} of \code{data.frame()} with the statistical results. This is
similar to \code{fetch_data("modeling_results")}.
}
\description{
This function is provided for convenience. It runs all the functions
required for computing the \code{modeling_results}. This can be useful for
finding marker genes on a new spatially-resolved transcriptomics dataset
and thus using it for \code{run_app()}. The results from this function can also be
used for performing spatial registration through \code{layer_stat_cor()} and
related functions of sc/snRNA-seq datasets.
}
\examples{
## Ensure reproducibility of example data
set.seed(20220907)

## Generate example data
sce <- scuttle::mockSCE()

## Add some sample IDs
sce$sample_id <- sample(LETTERS[1:5], ncol(sce), replace = TRUE)

## Add a sample-level covariate: age
ages <- rnorm(5, mean = 20, sd = 4)
names(ages) <- LETTERS[1:5]
sce$age <- ages[sce$sample_id]

## Add gene-level information
rowData(sce)$ensembl <- paste0("ENSG", seq_len(nrow(sce)))
rowData(sce)$gene_name <- paste0("gene", seq_len(nrow(sce)))

## Compute all modeling results
example_modeling_results <- registration_wrapper(
    sce,
    "Treatment", "sample_id", c("age"), "ensembl", "gene_name", "wrapper"
)
}
\seealso{
Other spatial registration and statistical modeling functions.: 
\code{\link{registration_block_cor}()},
\code{\link{registration_model}()},
\code{\link{registration_pseudobulk}()},
\code{\link{registration_stats_anova}()},
\code{\link{registration_stats_enrichment}()},
\code{\link{registration_stats_pairwise}()}
}
\concept{spatial registration and statistical modeling functions.}