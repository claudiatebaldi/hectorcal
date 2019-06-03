context('PCA functions')

library(magrittr)

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
    ncomp <- ncol(test_pca$rotation)

    p1 <- rep(0, ncomp); p1[1] <- 1        # expected result for runid==1
    p2 <- rep(0, ncomp); p2[1] <- 1; p2[2] <- 0.5 # expected result for runid==2
    p3 <- rep(NA, ncomp)               # runid==3 is real model output,
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
        ## data are equal to one another.  (Well, actually, this is only
        ## guaranteed if we kept all of the PCs for our test_pca data. In fact,
        ## we trimmed that set to avoid having a gigantic test set, so this is
        ## only approximately true.  As long as we keep enough components in the
        ## test set, the results will be sufficiently equal to satisfy the test.)
        ## However, we need to ensure that the
        ## input data is in canonical order.  We also need to filter out rows
        ## prior to the start year
        cd <- dplyr::arrange(climate_data, experiment, variable, year) %>%
            dplyr::filter(year >= min(test_pca$meta_data$histyear),
                          year <= max(test_pca$meta_data$year))
        expect_equal(recon_climate$value, cd$value,
                     info=paste("Reconstructed climate not equal to input for runid=",i))
    }

    ## Verify that the row-vector version of project_climate works
    cd <- test_climdata[test_climdata$runid==1, ]
    p <- project_climate(cd, test_pca)
    expect_equal(dim(p), c(1, ncomp))
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


test_that('summarize_pcdecomp works', {
    ## Create some mock-up data
    npc <- 50
    nmodel <- 10                   # changing this will change the test answers,
                                   # so don't
    pcs <-
        dplyr::bind_rows(
            lapply(1:nmodel,
                   function(modelid) {
                       pcnum <- 1:npc
                       value <- 1 + ((pcnum + modelid) %% nmodel) # number from 1-10
                       data.frame(model=as.character(modelid), PC=pcnum,
                                  value=value, stringsAsFactors=FALSE)
                   }))

    ## now the tests
    t1 <- summarize_pcdecomp(pcs)       # summarize full
    expect_equal(nrow(t1), npc)
    expect_equal(t1$mina, rep(1, npc))
    expect_equal(t1$maxb, rep(10, npc))
    expect_equal(t1$a10, rep(1.9, npc))
    expect_equal(t1$b90, rep(9.1, npc))
    expect_equal(t1$cmedian, rep(5.5, npc))
    expect_equal(t1$cmean, rep(5.5, npc))
    expect_equal(t1$PC, 1:npc)

    t2 <- summarize_pcdecomp(pcs, 10)   # First 10 only
    expect_equal(nrow(t2), 10)
    expect_equal(t2$PC, 1:10)
    expect_true(all(t1$mina==1))
    expect_true(all(t1$maxb==10))

    t3 <- summarize_pcdecomp(pcs, pcselect=c(2,4,6,8,10))
    expect_equal(nrow(t3), 5)
    expect_equal(t3$PC, c(2,4,6,8,10))

    t4 <- summarize_pcdecomp(pcs, 5, c(2,4,6,8,10)) # only 2 and 4
    expect_equal(nrow(t4), 2)
    expect_equal(t4$PC, c(2,4))
})
