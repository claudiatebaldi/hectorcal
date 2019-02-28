context('PCA functions')


test_that('Helper functions work', {

    test_df <- data.frame('test1' = rep(3, 3), 'test2' = rep(4, 3))

    # Test that the function works, nothing should be returned.
    expect_null(check_columns(input = test_df, req_cols = c('test1', 'test2')))

    # Test that the function throws an error
    expect_error(check_columns(input = test_df, req_cols = c('missing', 'test2')), "Missing columns: missing")

})

test_that('project_climate works', {

    # Load the test climate data and principal components.
    test_climdata <- readRDS('test_climdata.rds')
    test_pca <- readRDS('test_pca.rds')

    p1 <- rep(0, 50); p1[1] <- 1        # expected result for runid==1
    p2 <- rep(0, 50); p2[1] <- 1; p2[2] <- 0.5 # expected result for runid==2
    p3 <- rep(NA, 50)                          # runid==3 is real model output,
                                        # so we don't have an expected projection
                                        # for it.
    expected_proj <- list(p1, p2, p3)
    for(i in 1:3) {
        ## Pull out data from a single run to use as the input climate data.
        climate_data <- test_climdata[test_climdata$runid == i, ]

        ## Project the cliamte data on to the principal componet basis and then
        ## reconstruct the original climate data.
        projected_rslt <- project_climate(climate_data = climate_data,
                                          principal_components = test_pca, row_vector=FALSE)
        if(!any(is.na(expected_proj[[i]]))) {
            expect_equal(projected_rslt, expected_proj[[i]],
                         info=paste("Didn't get expected projection coefs for runid=",i))
        }
        recon_climate  <- reconstruct_climate(projected_rslt, principal_components = test_pca)

        ## We expect that the reconstructed climate data and the input climare
        ## data are equal to one another.  However, we need to ensure that the
        ## input data is in canonical order
        cd <- dplyr::arrange(climate_data, experiment, variable, year)
        expect_equal(recon_climate$value, cd$value,
                     info=paste("Reconstructed climate not equal to input for runid=",i),
                     tolerance=1e-3)    # The round-trip procedure can induce a
                                        # fair bit of roundoff error.
    }

    ## Verify that the row-vector version of project_climate works
    cd <- test_climdata[test_climdata$runid==1, ]
    p <- project_climate(cd, test_pca)
    expect_equal(dim(p), c(1, 50))
    expect_equal(as.vector(p), p1)
    rclim <- reconstruct_climate(p, test_pca)
    expect_equal(rclim$value, cd$value)

    ## Verify that setting ncomp in reconstruct_climate works
    c1 <- reconstruct_climate(p2, test_pca)
    c2 <- reconstruct_climate(p2, test_pca, 1)
    vdiff <- (c1$value - c2$value) / test_pca$scale
    expect_equivalent(vdiff, 0.5*test_pca$rotation[ ,2])

    test_PC <- test_pca    # avoid having to change names in all the tests below.
    # Make sure that the project climate function throws an error, when the climate input data and the principal components
    # do not have the same data structure.
    missing_var <- climate_data
    missing_var <- dplyr::filter(missing_var, variable == 'co2')
    expect_error(project_climate(missing_var, test_PC), 'nrow.cd. not equal to nrow.principal_components')

    missing_years <- climate_data
    missing_years$year <- c(climate_data$year + 10)
    expect_error(project_climate(missing_years, test_PC), 'nrow.cd. not equal to nrow.principal_components')

    missing_experiment            <- climate_data
    missing_experiment$experiment <- c('rcp85')
    expect_error(project_climate(missing_experiment, test_PC), 'experiment %in% climate_data')

    # pc missing columns
    missing_pc_rows <- test_PC
    missing_pc_rows$rotation <- missing_pc_rows$rotation[1:3, ]
    expect_error(project_climate(climate_data, missing_pc_rows), 'nrow.cd. not equal to nrow.principal_components')

})
