## Obtain the necessary data
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

## Obtain labels
annotate_registered_clusters(cor_stats_layer)


test_that("Stop with slash", {
  colnames(cor_stats_layer) <- gsub("ayer", "/", colnames(cor_stats_layer))
  expect_error(annotate_registered_clusters(cor_stats_layer))
})
