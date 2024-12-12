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

## add variable with one group
sce$batch <- "batch1"

## Add gene-level information
rowData(sce)$ensembl <- paste0("ENSG", seq_len(nrow(sce)))
rowData(sce)$gene_name <- paste0("gene", seq_len(nrow(sce)))


test_that(
    "warning for k=2 variable",
    example_modeling_results <- expect_warning(
        registration_wrapper(
            sce,
            var_registration = "Treatment",
            var_sample_id = "sample_id",
            covars = c("age"),
            gene_ensembl = "ensembl",
            gene_name = "gene_name",
            suffix = "wrapper"
        )
    )
)

#### Syntactic Variable Test ####

## catagorical var as int
sce$cluster_int <- sample(1:4, ncol(sce), replace = TRUE)
# sce$cluster_k <- paste0("k", sce$cluster_int)
sce$cluster_j <- paste0(sce$cluster_int,"j")
sce$cluster_l <- sample(c("L-1", "L2/3", "4L", "L5"), ncol(sce), replace = TRUE)

table(sce$cluster_j)

test_that("Numeric var_regisration throws warning",
    expect_warning(
        registration_wrapper(
            sce,
            var_registration = "cluster_int",
            var_sample_id = "sample_id",
            covars = c("age"),
            gene_ensembl = "ensembl",
            gene_name = "gene_name",
            suffix = "wrapper"
        )
    ))

test_that("Non-Syntactic throws warning",
    expect_warning(
        registration_wrapper(
            sce,
            var_registration = "cluster_l",
            var_sample_id = "sample_id",
            covars = c("age"),
            gene_ensembl = "ensembl",
            gene_name = "gene_name",
            suffix = "wrapper"
        )
    ))








