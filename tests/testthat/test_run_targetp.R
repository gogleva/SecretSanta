context("Check targetp")

test_that("targetp correctly responds to invalid inputs",
          {
            # prep input object
            aa <- readAAStringSet(system.file("extdata",
                                              "sample_prot_100.fasta",
                                              package = "SecretSanta"),
                                  use.names = TRUE)
            inp <- CBSResult(in_fasta = aa[1:10])
            
            # test with inp_object belonging to an incorrect class:
            expect_error(suppressMessages(targetp(input_obj = aa,
                                                  network = 'N',
                                                  run_mode = "starter")),
                         "input_obj does not belong to CBSResult superclass")
            
            # test starters with valid input options:
            
            expect_is(suppressMessages(targetp(inp,'N',run_mode = "starter")),
                      "TargetpResult")
            expect_is(suppressMessages(targetp(inp, 'P', run_mode = "starter")),
                      "TargetpResult")
            
            # piper:
            expect_error(suppressMessages(targetp(inp, 'P', run_mode = "piper")),
                      "out_fasta attribute is empty")
            
          
            # # test pipers with valid input options:
            
             s1_sp2 <- signalp(inp,
                               version = 2,
                               organism = 'euk',
                               run_mode = "starter")
             s2_tp <- targetp(s1_sp2,
                              network = 'P',
                              run_mode = "piper")
             expect_is(s2_tp, "TargetpResult")
             
             # test starter with empty in_fasta attribute
             
              emp <- SignalpResult()
              expect_error(targetp(emp,
                                   network = 'N',
                                   run_mode = "starter"),
                          "in_fasta attribute is empty")
             
            # # test piper with empty out_fasta attribute
            
              expect_error(targetp(emp,
                                   network = 'N',
                                   run_mode = "piper"),
                          "out_fasta attribute is empty")
            
            
            # # test invalid network type
            expect_error(suppressMessages(targetp(inp,
                                                  network = 'Plant',
                                                  run_mode = "starter")), 
                         "'arg' should be one of “P”, “N”")
            expect_error(suppressMessages(targetp(inp,
                                                  network = 'n',
                                                  run_mode = "starter")), 
                         "'arg' should be one of “P”, “N”")
            expect_error(suppressMessages(targetp(inp,
                                                  network = 'p',
                                                  run_mode = "starter")), 
                         "'arg' should be one of “P”, “N”")

        })
