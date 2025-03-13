
## Obtain the necessary data
## reference human pilot modeling results
if (!exists("modeling_results")) {
    modeling_results <- fetch_data(type = "modeling_results")
}

## query spatialDLPFC modeling results
query_modeling_results <- fetch_data(
    type = "spatialDLPFC_Visium_modeling_results"
)

## Compute the correlations
cor_stats_layer <- layer_stat_cor(
    stats = query_modeling_results$enrichment,
    modeling_results,
    model_type = "enrichment"
)

## Default plot with no annotations and defaults for ComplexHeatmap()
# default_plot <- layer_stat_cor_plot(cor_stats_layer)


## Add annotation
annotation_df <- annotate_registered_clusters(
    cor_stats_layer,
    confidence_threshold = .55
)

# anno_plot <- layer_stat_cor_plot(cor_stats_layer, annotation = annotation_df)

## test w/ numeric reference lables
cor_stats_numeric <- cor_stats_layer
colnames(cor_stats_numeric) <- c(13, 6:1)

annotation_df_numeric <- annotate_registered_clusters(
  cor_stats_numeric,
  confidence_threshold = .55
)

# layer_stat_cor_plot(cor_stats_numeric, annotation = annotation_df_numeric)


test_that("annotation checks work",{
  annotation_df$cluster <- gsub("Sp09", "",annotation_df$cluster)
  expect_error(create_annotation_matrix(annotation_df, cor_stats_layer))
})

test_that("annotation matrix works w/ short cluster names", {
  anno_matrix <- create_annotation_matrix(annotation_df, cor_stats_layer)
  anno_matrix_numeric <- create_annotation_matrix(annotation_df_numeric, cor_stats_numeric)
  colnames(anno_matrix_numeric) <- colnames(anno_matrix)
  
  expect_equal(anno_matrix, anno_matrix_numeric)
})

