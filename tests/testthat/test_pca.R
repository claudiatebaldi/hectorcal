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
    test_PCA <- get(load(list.files(system.file('tests/testthat', package = 'hectorcal'), 'PCA_climate_test.rda', full.names = TRUE)))
    test_PC  <- get(load(list.files(system.file('tests/testthat', package = 'hectorcal'), 'PC_climate_test.rda', full.names = TRUE)))

    # Pull out data from a single run to use as the input climate data.
    climate_data <- test_PCA[test_PCA$runid == 2, ]

    # Project the cliamte data on to the principal componet basis and then reconstruct the original climate data.
    projected_rslt <- project_climate(climate_data = climate_data, principal_components = test_PC)
    recon_climate  <- reconstruct_climate(projected_rslt, principal_components = test_PC)

    # We expect that the reconstructed climate data and the input climare data are equal to one another.
    expect_equal(recon_climate$value, climate_data$value)

    # Make sure that the project climate function throws an error, when the climate input data and the principal components
    # do not have the same data structure.
    missing_var <- climate_data
    missing_var$variable <- 'co2'
    expect_error(project_climate(missing_var, test_PC), 'climate_data and principal_components are required to have the same variables')

    missing_years <- climate_data
    missing_years$year <- c(climate_data$year + 10)
    expect_error(project_climate(missing_years, test_PC), 'climate_data and principal_components are required to have the same years')

    missing_experiment            <- climate_data
    missing_experiment$experiment <- c('rcp26', 'rcp85', 'rcp26', 'rcp26', 'rcp85', 'rcp26', 'rcp26', 'rcp85', 'rcp26', 'rcp45', 'rcp45')
    expect_error(project_climate(missing_experiment, test_PC), 'climate_data and principal_components are required to have the same experiments')

    # pc missing columns
    missing_pc_rows <- test_PC
    missing_pc_rows$rotation <- missing_pc_rows$rotation[1:3, ]
    expect_error(project_climate(climate_data, missing_pc_rows), 'There is a mismatch between the information in the climate_data and the principal_components arguments.')

})
