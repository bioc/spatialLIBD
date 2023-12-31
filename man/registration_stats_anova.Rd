% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/registration_stats_anova.R
\name{registration_stats_anova}
\alias{registration_stats_anova}
\title{Spatial registration: compute ANOVA statistics}
\usage{
registration_stats_anova(
  sce_pseudo,
  block_cor,
  covars = NULL,
  var_registration = "registration_variable",
  var_sample_id = "registration_sample_id",
  gene_ensembl = NULL,
  gene_name = NULL,
  suffix = ""
)
}
\arguments{
\item{sce_pseudo}{The output of \code{registration_pseudobulk()}.}

\item{block_cor}{A \code{numeric(1)} computed with \code{registration_block_cor()}.}

\item{covars}{A \code{character()} with names of sample-level covariates.}

\item{var_registration}{A \code{character(1)} specifying the \code{colData(sce_pseudo)}
variable of interest against which will be used for computing the relevant
statistics.}

\item{var_sample_id}{A \code{character(1)} specifying the \code{colData(sce_pseudo)}
variable with the sample ID.}

\item{gene_ensembl}{A \code{character(1)} specifying the \code{rowData(sce_pseudo)}
column with the ENSEMBL gene IDs. This will be used by \code{layer_stat_cor()}.}

\item{gene_name}{A \code{character(1)} specifying the \code{rowData(sce_pseudo)}
column with the gene names (symbols).}

\item{suffix}{A \code{character(1)} specifying the suffix to use for the
F-statistics column. This is particularly useful if you will run this
function more than once and want to be able to merge the results.}
}
\value{
A \code{data.frame()} with the ANOVA statistical results. This is
similar to \code{fetch_data("modeling_results")$anova}.
}
\description{
This function computes the gene ANOVA F-statistics (at least one group is
different from the rest). These F-statistics can be used for spatial
registration with \code{layer_stat_cor()} and related functions. Although, they
are more typically used for identifying ANOVA-marker genes.
}
\examples{
example("registration_block_cor", package = "spatialLIBD")
results_anova <- registration_stats_anova(sce_pseudo,
    block_cor, "age",
    gene_ensembl = "ensembl", gene_name = "gene_name", suffix = "example"
)
head(results_anova)

## Specifying `block_cor = NaN` then ignores the correlation structure
results_anova_nan <- registration_stats_anova(sce_pseudo,
    block_cor = NaN, "age",
    gene_ensembl = "ensembl", gene_name = "gene_name", suffix = "example"
)
head(results_anova_nan)

## Note that you can merge multiple of these data.frames if you run this
## function for different sets. For example, maybe you drop one group
## before pseudo-bulking if you know that there are many differences between
## that group and others. For example, we have dropped the white matter (WM)
## prior to computing ANOVA F-statistics.

## no covariates
results_anova_nocovar <- registration_stats_anova(sce_pseudo,
    block_cor,
    covars = NULL,
    gene_ensembl = "ensembl", gene_name = "gene_name", suffix = "nocovar"
)
head(results_anova_nocovar)

## Merge both results into a single data.frame, thanks to having different
## 'suffix' values.
results_anova_merged <- merge(results_anova, results_anova_nocovar)
head(results_anova_merged)
}
\seealso{
Other spatial registration and statistical modeling functions: 
\code{\link{registration_block_cor}()},
\code{\link{registration_model}()},
\code{\link{registration_pseudobulk}()},
\code{\link{registration_stats_enrichment}()},
\code{\link{registration_stats_pairwise}()},
\code{\link{registration_wrapper}()}
}
\concept{spatial registration and statistical modeling functions}
