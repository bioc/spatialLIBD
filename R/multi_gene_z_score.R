#' Combine multiple continuous variables by averaging Z scores
#'
#' To summarize multiple features, each is normalized to represent a Z-score.
#' Scores are averaged to return a single vector.
#'
#' @param cont_mat A \code{matrix()} with spots as rows and 2 or more continuous
#' variables as columns.
#'
#' @return A \code{numeric()} vector with one element per spot, summarizing the
#' multiple continuous variables.
#'
#' @author Nicholas J. Eagles
#' @import MatrixGenerics
#' @family functions for summarizing expression of multiple continuous variables simultaneously
#' @keywords internal
multi_gene_z_score <- function(cont_mat) {
    #   Z-score calculation requires at least 1 feature with nonzero variance.
    #   Verify this and drop any zero-variance features
    good_indices <- which(colSds(cont_mat, na.rm = TRUE) != 0)
    if (length(good_indices) < 1) {
        stop("After dropping features with no expression variation, no features were left. This error can occur when using data from only 1 spot.", call. = FALSE)
    }
    if (ncol(cont_mat) - length(good_indices) > 0) {
        warning(
            sprintf(
                "Dropping features(s) '%s' which have no expression variation",
                paste(colnames(cont_mat)[-good_indices], collapse = "', '")
            ),
            call. = FALSE
        )
    }
    cont_mat <- cont_mat[, good_indices, drop = FALSE]

    #   For each spot, average Z-scores across all features
    cont_z <- (cont_mat - rep(colMeans(cont_mat, na.rm = TRUE), each = nrow(cont_mat))) /
        rep(colSds(cont_mat, na.rm = TRUE), each = nrow(cont_mat))
    z_vec <- rowMeans(cont_z, na.rm = TRUE)

    return(z_vec)
}
