context("Check m-slicer")

test_that("m-slicer treats input objects correctly",
          {
            # prep input objects
            
            aa <- readAAStringSet(
              system.file("extdata",
                          "sample_prot_100.fasta",
                          package = "SecretSanta"),
              use.names = TRUE
            )
            inp <- CBSResult(in_fasta = aa[1:10])
            s1_sp2 <- signalp(inp,
                              version = 2,
                              organism = 'euk',
                              run_mode = "starter",
                              legacy_method = 'hmm')
            
            # test with valid inputs
            
            r1 <- m_slicer(aa, 100, run_mode = 'slice')
            expect_is(r1, 'AAStringSet')
            expect_true(length(r1) > length(aa))
            
            expect_error(
              m_slicer(aa,
                       100,
                       run_mode = 'rescue'),
              "Use run_mode 'slice' for an input object of AAStringSet class"
            )
            
            #test the rescue mode:
            r2 <- m_slicer(s1_sp2, 100, run_mode = 'rescue')
            expect_is(r2, 'AAStringSet')
            
          })
