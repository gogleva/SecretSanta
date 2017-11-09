context("Check external dependencies")

test_that("manage_paths returns correct object and messages",
          {
            s_file <- system.file("extdata", "sample_paths",
                                  package = "SecretSanta")
            
            # check with in_path = TRUE (assume that all dependecies should be
            # accessible via $PATH variable)
            
            # this test is not crucial, can be suspended
            expect_message(
              manage_paths(
                in_path = TRUE,
                test_tool = 'all',
                'checking dependencies acessible via $PATH'
              )
            )
            
            # this test is crucial
            r1 <- manage_paths(in_path = TRUE, test_tool = 'all')
            expect_is(r1, 'list')
            expect_true(r1["tests"] == TRUE)
            expect_true(r1["in_path"] == TRUE)
            expect_true(is.na(r1["path_tibble"]))
            
            # check with in_path = FALSE and provide in_path file
            
            # this test is not crucial, can be suspended
            expect_message(
              manage_paths(
                in_path = FALSE,
               test_tool = 'all',
                path_file = s_file
              ),
              'All paths are valid'
            )
            
            # this test is crucial:
            r2 <- manage_paths(in_path = FALSE,
                              test_tool = 'all',
                               path_file = s_file)
            
            expect_is(r2, 'list')
            expect_true(r2["tests"] == TRUE)
            expect_true(r2["in_path"] == FALSE)
            expect_is(r2$path_tibble, 'tbl')
            expect_equal(ncol(r2$path_tibble), 3)
          })
