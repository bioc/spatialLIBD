## Obtain the necessary data
if (!exists("modeling_results")) {
    modeling_results <- fetch_data(type = "modeling_results")
}
if (!exists("sce_layer")) sce_layer <- fetch_data(type = "sce_layer")

## enrichment top 10 genes
enrich_top10 <- sig_genes_extract(
    n = 10,
    modeling_results = modeling_results,
    model_type = "enrichment",
    sce_layer = sce_layer
)

test_that("Extract expected number of enrichment genes", {
    n_tests <- sum(grepl("t_stat", colnames(modeling_results$enrichment)))

    expect_equal(nrow(enrich_top10), n_tests * 10)
})

test_that("gene_name is flexible", {
    ## rename gene_name column
    rowData_names <- colnames(rowData(sce_layer))
    colnames(rowData(sce_layer))[grep("gene_name", rowData_names)] <- "Symbol"

    ## Enrichment top 10 genes - use "Symbol" as gene_name
    enrich_top10_change_name <- sig_genes_extract(
        n = 10,
        model = "enrichment",
        modeling_results = modeling_results,
        sce_layer = sce_layer,
        gene_name = "Symbol"
    )

    expect_equal(enrich_top10, enrich_top10_change_name)
})
