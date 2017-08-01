context("Check paths")

test_that("manage_paths function works as expected",
          {
          expect_message(manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta")),
                           'All paths are valid')
          
          result <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
          expect_is(result,
                    c('tbl_df', 'tbl'))
          
          })