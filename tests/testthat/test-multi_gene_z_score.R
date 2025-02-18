test_that(
    "multi_gene_z_score",
    {
        #   With two good columns but 1 zero-variance column, the zero-variance
        #   column should be dropped with a warning
        cont_mat <- matrix(c(1, 0, 3, 3, 2, -5), ncol = 3)
        colnames(cont_mat) <- c("good1", "bad", "good2")
        expect_warning(
            multi_gene_z_score(cont_mat),
            "Dropping features\\(s\\) 'bad' which have no expression variation"
        )

        #   NAs should be correctly removed from columns (as long as 2 non-NAs remain
        #   in at least 1 column), and the result should have no NAs
        cont_mat <- matrix(c(1, NA, 3, NA, 2, 0), ncol = 2)
        colnames(cont_mat) <- c("good1", "good2")
        expect_equal(any(is.na(multi_gene_z_score(cont_mat))), FALSE)

        #   With only one good column, the result should simply be the
        #   Z-score-normalized good column. A warning should indicate which
        #   columns were dropped
        cont_mat <- matrix(c(1, NA, 3, 4, 2, 2), ncol = 3)
        colnames(cont_mat) <- c("bad1", "good", "bad2")

        temp <- c(3, 4)
        expected_result <- (temp - mean(temp)) / sd(temp)

        expect_warning(
            {
                actual_result <- multi_gene_z_score(cont_mat)
            },
            "Dropping features\\(s\\) 'bad1', 'bad2' which have no expression variation"
        )
        expect_equal(actual_result, expected_result)

        #   An error should be thrown if no columns have variation
        cont_mat <- matrix(c(1, 1, 0, 0, 2, 2), ncol = 3)
        colnames(cont_mat) <- c("bad1", "bad2", "bad3")
        expect_error(
            multi_gene_z_score(cont_mat),
            "^After dropping features with no expression variation"
        )

        #   Now actually verify the Z-score-based calculation is correct
        cont_mat = matrix(c(1, 2, 3, -3, 0, 3, 7.5, 9, 10.5), nrow = 3)
        expect_equal(multi_gene_z_score(cont_mat), c(-1, 0, 1))
    }
)
