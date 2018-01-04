context("Check signalp")

test_that("signalp correctly responds to invalid inputs",
          {
          # prep input object
          aa <- readAAStringSet(system.file("extdata", 
                                            "sample_prot_100.fasta",
                                            package = "SecretSanta"),
                                use.names = TRUE)
          inp <- SignalpResult(in_fasta = aa[1:10])
          
          # test with dna set in the in_fasta slot:
          dna <- readAAStringSet(system.file("extdata",
                                             "sample_dna.fasta",
                                             package = "SecretSanta"))
          expect_error(CBSResult(in_fasta = dna)) 
          
          # test with inp_object belonging to an incorrect class:
          expect_error(suppressMessages(signalp(aa,
                                                organism = 'euk',
                                                version = 2, 
                                                run_mode = "starter")),
                       "input_object does not belong to CBSResult superclass")
          expect_error(suppressMessages(signalp(inp,
                                             version = 2,
                                             organism = 'euk',
                                             run_mode = "starter")),
                    'missing argument: legacy_method')
          
          # test starters with valid input options:
          
          expect_is(suppressMessages(signalp(inp,
                                             version = 2,
                                             organism = 'euk',
                                             run_mode = "starter",
                                             legacy_method = 'hmm')),
                                             "SignalpResult")
          expect_is(suppressMessages(signalp(inp,
                                             version = 3,
                                             organism = 'euk',
                                             run_mode = "starter",
                                             legacy_method = 'hmm')),
                                             "SignalpResult")
          
          expect_is(suppressMessages(signalp(inp,
                                             version = 4,
                                             organism = 'euk',
                                             run_mode = "starter")),
                                             "SignalpResult")
          
          expect_is(suppressMessages(signalp(inp,
                                             version = 3.0,
                                             organism = 'euk',
                                             run_mode = "starter",
                                             legacy_method = 'hmm')), 
                                             "SignalpResult")
          
          # test pipers with valid input options:
          
          s1_sp2 <- signalp(inp,
                            version = 2, 
                            organism = 'euk',
                            run_mode = "starter",
                            legacy_method = 'hmm')
          expect_is(s1_sp2, "SignalpResult")
          
          s2_sp3 <- signalp(s1_sp2,
                            version = 3,
                            organism = 'euk',
                            run_mode = "piper",
                            legacy_method = 'hmm')
          expect_is(s2_sp3, "SignalpResult")
          s3_sp4 <- signalp(s2_sp3,
                            version = 4,
                            organism = 'euk',
                            run_mode = "piper")
          expect_is(s3_sp4, "SignalpResult")
          
          s4_sp41 <- signalp(s2_sp3,
                            version = 4.1,
                            organism = 'euk',
                            run_mode = "piper")
          expect_is(s4_sp41, "SignalpResult")
          
          # test starter with empty in_fasta attribute
          
          emp <- SignalpResult()
          expect_error(signalp(emp,
                               version = 2,
                               organism = "euk",
                               run_mode = "starter"),
                       "in_fasta attribute is empty")
          
          # test piper with empty out_fasta attribute
          expect_error(signalp(emp,
                               version = 2,
                               organism = "euk",
                               run_mode = "piper",
                               paths = my_pa),
                       "out_fasta attribute is empty")
          
          expect_error(signalp(inp,
                               version = 2,
                               organism = "euk",
                               run_mode = "piper"),
                       "out_fasta attribute is empty")
          
          # test invalid versions
          expect_error(suppressMessages(signalp(
            inp,
            version = 5, 
            organism = 'euk',
            run_mode = "starter")))
            
          expect_error(suppressMessages(signalp(
            inp,
            version = 1,
            organism = 'euk',
            run_mode = "starter",
            paths = my_pa))) 
 
          # test invalid organism
          expect_error(suppressMessages(signalp(inp,
                                                version = 3,
                                                organism = 'bacteria',
                                                run_mode = "starter")))
          expect_error(suppressMessages(signalp(inp,
                                                version = 3, 
                                                organism = 'gram',
                                                run_mode = "starter"))) 
          
          
          ## test sensitive mode when signalp4.1 is used:
          
          
          expect_is(signalp(inp,
                  version = 3,
                  organism = 'euk',
                  run_mode = 'starter',
                  sensitive = TRUE,
                  legacy_method = 'hmm'), 
                  "SignalpResult")
          
          expect_message(signalp(inp,
                  version = 4,
                  organism = 'euk',
                  run_mode = 'starter',
                  sensitive = TRUE,
                  legacy_method = 'hmm'),
                  'version > 3, legacy_method is not required')
          
          expect_equal(nrow(getSPtibble(signalp(inp = CBSResult(in_fasta = aa),
                  version = 4,
                  organism = 'euk',
                  run_mode = 'starter',
                  sensitive = TRUE))), 7) 
          
          expect_equal(nrow(getSPtibble(signalp(inp = CBSResult(in_fasta = aa),
                  version = 4,
                  organism = 'euk',
                  run_mode = 'starter',
                  sensitive = FALSE))), 5)

        })
