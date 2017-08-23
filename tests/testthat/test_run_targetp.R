context("Check targetp")

test_that("targetp correctly responds to invalid inputs",
          {
            # prep input object
            my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
            inp <- SignalpResult()
            aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"),
                                  use.names = TRUE)
            inp <- setInfasta(inp, aa[1:10])
            
            # test with inp_object belonging to an incorrect class:
            expect_error(suppressMessages(targetp(input_object = aa, network_type = 'N',
                                                  run_mode = "starter", paths = my_pa)),
                         "input_object does not belong to CBSResult superclass")
            
            # test starters with valid input options:
            
            expect_is(suppressMessages(targetp(inp, 'N', run_mode = "starter", paths = my_pa)),
                      "TargetpResult")
            expect_is(suppressMessages(targetp(inp, network_type = 'P', run_mode = "starter", paths = my_pa)),
                      "TargetpResult")
            
            # piper:
            expect_error(suppressMessages(targetp(inp, network_type = 'P', run_mode = "piper", paths = my_pa)),
                      "out_fasta attribute is empty")
            
          
            # # test pipers with valid input options:
            
             s1_sp2 <- signalp(inp, version = 2, 'euk', run_mode = "starter", paths = my_pa)
             s2_tp <- targetp(s1_sp2, network_type = 'P', run_mode = "piper", paths = my_pa)
             expect_is(s2_tp, "TargetpResult")
             
            # # test starter with empty in_fasta attribute
             
              emp <- SignalpResult()
              expect_error(targetp(emp, network_type = 'N', run_mode = "starter", paths = my_pa),
                          "in_fasta attribute is empty")
             
            # # test piper with empty out_fasta attribute
            
              expect_error(targetp(emp, network_type = 'N', run_mode = "piper", paths = my_pa),
                          "out_fasta attribute is empty")
            
            
            # # test invalid network type
            expect_error(suppressMessages(targetp(inp, network_type = 'Plant', run_mode = "starter", paths = my_pa)), 
                          "Specified network_type is invalid.")
            expect_error(suppressMessages(targetp(inp, network_type = 'n', run_mode = "starter", paths = my_pa)), 
                         "Specified network_type is invalid.")
            expect_error(suppressMessages(targetp(inp, network_type = 'p', run_mode = "starter", paths = my_pa)), 
                         "Specified network_type is invalid.")

            # # test invalid run mode
            expect_error(suppressMessages(targetp(inp, network_type = N, run_mode = "start", paths = my_pa)), 
                          "Run mode is invalid. Please use 'starter' to initiate prediction pipelie or 'piper' to continue")
            expect_error(suppressMessages(targetp(inp, network_type = N, run_mode = "pipe", paths = my_pa)), 
                         "Run mode is invalid. Please use 'starter' to initiate prediction pipelie or 'piper' to continue")
          })
