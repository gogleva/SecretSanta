context("Check signalp")

test_that("signalp correctly responds to invalid inputs",
          {
          # prep input object
          my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
          inp <- SignalpResult()
          aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"),
                                use.names = TRUE)
          inp <- setInfasta(inp, aa[1:10])
          
          # test with dna set in the in_fasta slot:
          dna <- readAAStringSet(system.file("extdata", "sample_dna.fasta", package = "SecretSanta"))
          expect_error(CBSResult(in_fasta = dna)) 
          
          # test with inp_object belonging to an incorrect class:
          expect_error(suppressMessages(signalp(aa, version = 2, 'euk', run_mode = "starter", paths = my_pa)),
                    "input_object does not belong to CBSResult superclass")
          
          # test starters with valid input options:
          
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
          
          # test pipers with valid input options:
          
          s1_sp2 <- signalp(inp, version = 2, 'euk', run_mode = "starter", paths = my_pa)
          expect_is(s1_sp2, "SignalpResult")
          s2_sp3 <- signalp(s1_sp2, version = 3, 'euk', run_mode = "piper", paths = my_pa)
          expect_is(s2_sp3, "SignalpResult")
          s3_sp4 <- signalp(s2_sp3, version = 4, 'euk', run_mode = "piper", paths = my_pa)
          expect_is(s3_sp4, "SignalpResult")
          
          # test starter with empty in_fasta attribute
          
          emp <- SignalpResult()
          expect_error(signalp(emp, version = 2, "euk", run_mode = "starter", paths = my_pa),
                       "in_fasta attribute is empty")
          
          # test piper with empty out_fasta attribute
          expect_error(signalp(emp, version = 2, "euk", run_mode = "piper", paths = my_pa),
                       "out_fasta attribute is empty")
          
          expect_error(signalp(inp, version = 2, "euk", run_mode = "piper", paths = my_pa),
                       "out_fasta attribute is empty")
          
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
          
          # test invalid run mode
          expect_error(suppressMessages(signalp(inp, version = 3, 'gram-', run_mode = "start", paths = my_pa)), 
                       "Run mode is invalid. Please use 'starter' to initiate prediction pipelie or 'piper' to continue")
          })
