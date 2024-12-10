

#' Visualize the correlation of layer modeling t-statistics with ComplexHeatmap
#' @param query_colors named vector of colors for query row annotations
#' @param reference_colors named vector of colors for reference column annotations
#'
#' @inheritParams layer_stat_cor_plot
#'
#' @return ComplexHeatmap plot of t-stat correlations
#' @export
#'
#' @examples
#' ## Obtain the necessary data
#' ## reference human pilot modeling results
#' if (!exists("modeling_results")) {
#'     modeling_results <- fetch_data(type = "modeling_results")
#' }
#' 
#' ## querey spatailDLPFC modeling 
#' query_modeling_results <- fetch_data(type = "spatialDLPFC_Visium_modeling_results")
#' 
#' ## extract t-statics and rename
#' registration_t_stats <- query_modeling_results$enrichment[, grep("^t_stat", colnames(query_modeling_results$enrichment))]
#' colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))
#'
#' ## Compute the correlations
#' cor_stats_layer <- layer_stat_cor(
#'     stats = registration_t_stats,
#'     modeling_results,
#'     model_type = "enrichment"
#' )
#'
#' ## Visualize the correlation matrix
#' 
#' ## most basic
#' layer_stat_cor_plot_complex(cor_stats_layer)
#' 
#' ## add colors
#' ## add libd_layer_colors to refrence Human Pilot layers
#' layer_stat_cor_plot_complex(cor_stats_layer, reference_colors = libd_layer_colors)
#' 
#' ## supply polychrome colors to query clusters
#' cluster_colors <-  Polychrome::palette36.colors(nrow(cor_stats_layer))
#' names(cluster_colors) <- rownames(cor_stats_layer)
#'
#' layer_stat_cor_plot_complex(cor_stats_layer, 
#'                             query_colors = cluster_colors,
#'                             reference_colors = libd_layer_colors)
#' 
#' ## Apply additional ComplexHeatmap param
#' layer_stat_cor_plot_complex(cor_stats_layer, cluster_rows = FALSE)
#' 
layer_stat_cor_plot_complex <- function(cor_stats_layer,
                                        theSeq = seq(min(cor_stats_layer), max(cor_stats_layer), by = 0.01),
                                        my.col = grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "PRGn"))(length(theSeq)),
                                        query_colors = NULL,
                                        reference_colors = NULL,
                                        ...
                                        ){
  
  # ## query annotations on row
  if(!is.null(query_colors)){
    
    stopifnot(all(rownames(cor_stats_layer) %in% names(query_colors)))
    query_colors <- query_colors[rownames(cor_stats_layer)]
    
    
    query_row_annotation <- ComplexHeatmap::rowAnnotation(
      " " = rownames(cor_stats_layer),
      col = list(" " = query_colors),
      show_legend = FALSE)
    
  } else query_row_annotation <- NULL
  
  ## reference annotation on bottom
  if(!is.null(reference_colors)){
    stopifnot(all(colnames(cor_stats_layer) %in% names(reference_colors)))
    reference_colors <- reference_colors[colnames(cor_stats_layer)]
    
    ref_col_annotation <-  ComplexHeatmap::columnAnnotation(
      " " = colnames(cor_stats_layer),
      col = list(" " = reference_colors),
      show_legend = FALSE
    )
  } else ref_col_annotation <- NULL
  
  
  
  ## plot heatmap
  ComplexHeatmap::Heatmap(
    matrix = cor_stats_layer,
    col = my.col,
    bottom_annotation = ref_col_annotation,
    right_annotation = query_row_annotation,
    ...
  )
  
  
}
