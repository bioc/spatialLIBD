#' Visualize the correlation of layer modeling t-statistics with ComplexHeatmap
#'
#' This function makes a ComplexHeatmap from the correlation matrix
#' between a reference and query modeling statistics from [layer_stat_cor()].
#' For example, between the query statistics from a set of cell cluster/types
#' derived from scRNA-seq or snRNA-seq data (among other types) and the
#' reference layer statistics from the Human DLPFC Visium data (when using the
#' default arguments).
#'
#' Includes functionality to add color annotations,
#' (helpful to match to colors in Visium spot plots), and annotations from
#' [annotate_registered_clusters()].
#'
#' @param cor_stats_layer The output of [layer_stat_cor()].
#' @param color_max A `numeric(1)` specifying the highest correlation value for
#' the color scale (should be between 0 and 1).
#' @param color_min A `numeric(1)` specifying the lowest correlation value for
#' the color scale (should be between 0 and -1).
#' @param color_scale A `character(3)` vector specifying the color scale for the
#' fill of the heatmap. The first value is used for `color_min`, the second one
#' for zero, and the third for `color_max`.
#' @param query_colors named `character` vector of colors, Adds colors to query
#' row annotations.
#' @param reference_colors named `character` vector of colors, Adds colors to
#' reference column annotations.
#' @param annotation annotation data.frame output of
#' [annotate_registered_clusters()], adds 'X' for good confidence annotations,
#' '*' for poor confidence.
#' @param ... Additional parameters passed to
#' [ComplexHeatmap::Heatmap()][ComplexHeatmap::Heatmap()] such as `cluster_rows`
#' and `cluster_columns`.
#'
#'
#' @return ([Heatmap-class][ComplexHeatmap::Heatmap-class]) plot of t-stat
#' correlations
#' @export
#' @author Louise Huuki-Myers
#' @family Layer correlation functions
#'
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap columnAnnotation rowAnnotation Heatmap
#'
#' @examples
#' ## Obtain the necessary data
#' ## reference human pilot modeling results
#' if (!exists("modeling_results")) {
#'     modeling_results <- fetch_data(type = "modeling_results")
#' }
#'
#' ## query spatialDLPFC modeling results
#' query_modeling_results <- fetch_data(
#'     type = "spatialDLPFC_Visium_modeling_results"
#' )
#'
#' ## Compute the correlations
#' cor_stats_layer <- layer_stat_cor(
#'     stats = query_modeling_results$enrichment,
#'     modeling_results,
#'     model_type = "enrichment"
#' )
#'
#' ## Visualize the correlation matrix
#'
#' ## Default plot with no annotations and defaults for ComplexHeatmap()
#' layer_stat_cor_plot(cor_stats_layer)
#'
#' ## add colors
#' ## add libd_layer_colors to reference Human Pilot layers
#' layer_stat_cor_plot(cor_stats_layer, reference_colors = libd_layer_colors)
#'
#' ## obtain colors for the query clusters
#' cluster_colors <- get_colors(clusters = rownames(cor_stats_layer))
#' layer_stat_cor_plot(cor_stats_layer,
#'     query_colors = cluster_colors,
#'     reference_colors = libd_layer_colors
#' )
#'
#' ## Apply additional ComplexHeatmap param
#' layer_stat_cor_plot(cor_stats_layer,
#'     cluster_rows = FALSE,
#'     cluster_columns = FALSE
#' )
#'
#' ## Add annotation
#' annotation_df <- annotate_registered_clusters(
#'     cor_stats_layer,
#'     confidence_threshold = .55
#' )
#' layer_stat_cor_plot(cor_stats_layer, annotation = annotation_df)
#'
#' ## All together
#' layer_stat_cor_plot(
#'     cor_stats_layer,
#'     query_colors = cluster_colors,
#'     reference_colors = libd_layer_colors,
#'     annotation = annotation_df,
#'     cluster_rows = FALSE,
#'     cluster_columns = FALSE
#' )
#'
layer_stat_cor_plot <- function(cor_stats_layer,
    color_max = max(cor_stats_layer),
    color_min = min(cor_stats_layer),
    color_scale = c("#762A83", "#F7F7F7", "#1B7837"),
    query_colors = NULL,
    reference_colors = NULL,
    annotation = NULL,
    ...) {
    ## define color pallet
    stopifnot(color_min < color_max)
    stopifnot(color_min < 0)
    stopifnot(length(color_scale) == 3)
    my.col <- circlize::colorRamp2(c(color_min, 0, color_max), color_scale)

    # ## query annotations on row
    if (!is.null(query_colors)) {
        stopifnot(all(rownames(cor_stats_layer) %in% names(query_colors)))
        query_colors <- query_colors[rownames(cor_stats_layer)]


        query_row_annotation <- ComplexHeatmap::rowAnnotation(
            " " = rownames(cor_stats_layer),
            col = list(" " = query_colors),
            show_legend = FALSE
        )
    } else {
        query_row_annotation <- NULL
    }

    ## reference annotation on bottom
    if (!is.null(reference_colors)) {
        stopifnot(all(colnames(cor_stats_layer) %in% names(reference_colors)))
        reference_colors <- reference_colors[colnames(cor_stats_layer)]

        ref_col_annotation <- ComplexHeatmap::columnAnnotation(
            " " = colnames(cor_stats_layer),
            col = list(" " = reference_colors),
            show_legend = FALSE
        )
    } else {
        ref_col_annotation <- NULL
    }

    ## add annotation
    if (!is.null(annotation)) {
        anno_matrix <- create_annotation_matrix(annotation, cor_stats_layer)

        ## plot heatmap
        return(
            ComplexHeatmap::Heatmap(
                matrix = cor_stats_layer,
                col = my.col,
                name = "Cor",
                bottom_annotation = ref_col_annotation,
                right_annotation = query_row_annotation,
                cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(anno_matrix[i, j], x, y, gp = gpar(fontsize = 10))
                },
                ...
            )
        )
    }

    ## plot heatmap
    return(
        ComplexHeatmap::Heatmap(
            matrix = cor_stats_layer,
            col = my.col,
            name = "Cor",
            bottom_annotation = ref_col_annotation,
            right_annotation = query_row_annotation,
            ...
        )
    )
}

create_annotation_matrix <- function(annotation_df, cor_stats_layer) {
    anno_list <- lapply(
        rownames(cor_stats_layer),
        function(cluster) {
            # look up confidence
            confidence <- annotation_df[match(cluster, annotation_df$cluster), "layer_confidence"]
            sym <- ifelse(confidence == "good", "X", "*")
            # match annotations
            anno <- annotation_df[match(cluster, annotation_df$cluster), "layer_label"]
            return(ifelse(unlist(lapply(colnames(cor_stats_layer), grepl, anno)), sym, ""))
        }
    )

    anno_matrix <- t(data.frame(anno_list))
    rownames(anno_matrix) <- rownames(cor_stats_layer)
    colnames(anno_matrix) <- colnames(cor_stats_layer)

    return(anno_matrix)
}
