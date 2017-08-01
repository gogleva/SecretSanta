context("Check external dependencies")

test_that("manage_paths returns correct object and messages",
          {
          expect_message(manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta")),
                         'All paths are valid')
          expect_message(manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta")),
                           'signalp2 test run completed')
          expect_message(manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta")),
                         'signalp3 test run completed')
          expect_message(manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta")),
                         'signalp4 test run completed')
          expect_message(manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta")),
                         'targetp test run completed')
          expect_message(manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta")),
                         'tmhmm test run completed')
          expect_message(manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta")),
                         'wolfpsort test run completed')
          result <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
          expect_is(result, c('tbl_df', 'tbl'))
          expect_equal(ncol(result), 3)
          })