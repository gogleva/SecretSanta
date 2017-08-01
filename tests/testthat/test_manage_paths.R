context("Check external dependencies")

test_that("manage_paths returns correct object and messages",
          {
          s_file <- system.file("extdata", "sample_paths", package = "SecretSanta")
          
          expect_message(manage_paths(s_file), 'All paths are valid')
          expect_message(manage_paths(s_file), 'signalp2 test run completed')
          expect_message(manage_paths(s_file), 'signalp3 test run completed')
          expect_message(manage_paths(s_file), 'signalp4 test run completed')
          expect_message(manage_paths(s_file), 'targetp test run completed')
          expect_message(manage_paths(s_file), 'tmhmm test run completed')
          expect_message(manage_paths(s_file), 'wolfpsort test run completed')
          
          result <- manage_paths(s_file)
          expect_is(result, c('tbl_df', 'tbl'))
          expect_equal(ncol(result), 3)
          })