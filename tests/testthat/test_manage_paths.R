context("Check paths")

test_that("manage_paths fucntion works as expected",
          {
          expect_message(manage_paths(system.file("extdata", "sample_paths", package="SecretSanta")),
                       'All paths are valid')  
          })