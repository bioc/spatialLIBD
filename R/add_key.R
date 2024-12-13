#' Create a unique spot identifier
#'
#' This function adds `spe$key` to a
#' [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class] object
#' which is unique across all spots.
#'
#' @inheritParams add10xVisiumAnalysis
#' @param overwrite A `logical(1)` indicating whether to overwrite the `spe$key`.
#'
#' @return A
#' [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class] object
#' with `key` added to the `colData(spe)` that is unique across all spots.
#' @export
#'
#' @examples
#' if (enough_ram()) {
#'     ## Obtain the necessary data
#'     if (!exists("spe")) spe <- fetch_data("spe")
#'
#'     ## This object already has a 'key'
#'     head(spe$key)
#'
#'     ## We can clean it
#'     spe$key_original <- spe$key
#'     spe$key <- NULL
#'
#'     ## and then add it back
#'     spe <- add_key(spe)
#'     head(spe$key)
#'
#'     ## Note that the original 'key' order was 'sample_id'_'barcode' and we'
#'     ## have since changed it to 'barcode'_'sample_id'.
#'
#'     ## Below we restore the original 'key'
#'     spe$key <- spe$key_original
#'     spe$key_original <- NULL
#'     head(spe$key)
#' }
add_key <- function(spe, overwrite = TRUE) {
    if ("key" %in% colnames(colData(spe))) {
        if (overwrite) {
            message(
                "Overwriting 'spe$key'. Set 'overwrite = FALSE' if you do not want to overwrite it."
            )
            spe$key <- paste0(colnames(spe), "_", spe$sample_id)
            stopifnot(!any(duplicated(spe$key)))
        } else if (any(duplicated(spe$key))) {
            warning(
                "'spe$key' already exists and is not unique. Set 'overwrite = TRUE' to replace 'spe$key' with unique values.",
                call. = FALSE
            )
        }
    } else {
        spe$key <- paste0(colnames(spe), "_", spe$sample_id)
        stopifnot(!any(duplicated(spe$key)))
    }
    return(spe)
}
