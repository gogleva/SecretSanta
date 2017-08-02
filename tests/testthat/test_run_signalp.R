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
          starter_sp2 <- signalp(inp, version = 2, 'euk', run_mode = "starter", paths = my_pa)
          expect_is(starter_sp2, "SignalpResult")
          starter_sp3 <- signalp(inp, version = 3, 'euk', run_mode = "starter", paths = my_pa)
          expect_is(starter_sp3, "SignalpResult")
          starter_sp4 <- signalp(inp, version = 4, 'euk', run_mode = "starter", paths = my_pa)
          expect_is(starter_sp4, "SignalpResult")
          })

test_that("signalp outputs a valid result object",
          {
            
          })