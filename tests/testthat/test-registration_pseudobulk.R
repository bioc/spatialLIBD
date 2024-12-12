example("registration_pseudobulk", package = "spatialLIBD")

sce$batches <- sample(1:3, ncol(sce), replace = TRUE)
test_that("NA check works", {
    expect_error(
        registration_pseudobulk(sce, "Treatment", "sample_id", c("age", "batches")),
        "has all NAs after pseudo-bulking"
    )
    expect_error(
        registration_pseudobulk(sce, "CellCyle", "sample_id", c("age", "Treatment")),
        "var_registration"
    )
})


#### Syntactic Variable Test ####
set.seed(20220907) ## Ensure reproducibility of example data
sce <- scuttle::mockSCE()
## Add some sample IDs
sce$sample_id <- sample(LETTERS[1:5], ncol(sce), replace = TRUE)

## Add a sample-level covariate: age
ages <- rnorm(5, mean = 20, sd = 4)
names(ages) <- LETTERS[1:5]
sce$age <- ages[sce$sample_id]

## add variable with one group
sce$batch <- "batch1"

## non-syntactic inputs
sce$cluster_int <- sample(1:4, ncol(sce), replace = TRUE)
# sce$cluster_k <- paste0("k", sce$cluster_int)
sce$cluster_j <- paste0(sce$cluster_int,"j")
sce$cluster_l <- sample(c("L-1", "L2/3", "4L", "L5"), ncol(sce), replace = TRUE)

test_that("warn for numeric var_registration",
          expect_warning(registration_pseudobulk(sce,
                                  var_registration = "cluster_int",
                                  var_sample_id = "sample_id", 
                                  covars = c("age"), 
                                  min_ncells = NULL))
)


test_that("warn for non-syntactic var_registration",
          expect_warning(registration_pseudobulk(sce,
                                                 var_registration = "cluster_l",
                                                 var_sample_id = "sample_id", 
                                                 covars = c("age"), 
                                                 min_ncells = NULL))
)

