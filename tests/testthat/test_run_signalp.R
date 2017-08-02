context("Check signalp")

test_that("signalp correctly responds to invalid inputs",
          {
          # prep input object
          my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
          inp <- SignalpResult()
          aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"),
                                use.names = TRUE)
          inp <- setInfasta(inp, aa)
          
          # test with valid input options:
          
          expect_is(suppressMessages(signalp(inp, version = 2, 'euk', run_mode = "starter", paths = my_pa)),
                    "SignalpResult")
          expect_is(suppressMessages(signalp(inp, version = 3, 'euk', run_mode = "starter", paths = my_pa)),
                    "SignalpResult")
          expect_is(suppressMessages(signalp(inp, version = 4, 'euk', run_mode = "starter", paths = my_pa)),
                  "SignalpResult")
          expect_is(suppressMessages(signalp(inp, version = 3, 'Euk', run_mode = "starter", paths = my_pa)), 
                    "SignalpResult")
          expect_is(suppressMessages(signalp(inp, version = 3.0, 'Euk', run_mode = "starter", paths = my_pa)), 
                    "SignalpResult")
          
          
          # test invalid versions
          expect_error(suppressMessages(signalp(inp, version = 5, 'euk', run_mode = "starter", paths = my_pa)), 
                       'Input signalp version or specified organism type are invalid.')
          expect_error(suppressMessages(signalp(inp, version = 1, 'euk', run_mode = "starter", paths = my_pa)), 
                       'Input signalp version or specified organism type are invalid.')
          
          # test invalid organism
          expect_error(suppressMessages(signalp(inp, version = 3, 'bacteria', run_mode = "starter", paths = my_pa)), 
                       "Input signalp version or specified organism type are invalid.")
          expect_error(suppressMessages(signalp(inp, version = 3, 'gram', run_mode = "starter", paths = my_pa)), 
                       "Input signalp version or specified organism type are invalid.")
          
          
          })

test_that("signalp outputs a valid result object",
          {
            
          })