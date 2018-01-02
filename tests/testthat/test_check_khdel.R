context("Check ER-retention motif detection")

test_that("terminal KHDEL/HDEL motifs are detected",
          {
            # prep inputs:
            aa <- readAAStringSet(
              system.file("extdata",
                          "sample_prot_100.fasta",
                          package = "SecretSanta"),
              use.names = TRUE
            )
            inp <- CBSResult(in_fasta = aa[1:10])
            
            # check valid inputs, starter mode on a CBSResult object:

            # check object with empty out_fasta:
            expect_error(check_khdel(inp, pattern = 'strict'),
              'the input object contains empty out_fasta slot'
            )

            # chek outputs of signalp
            sp <-
              signalp(inp,
                      version = 2,
                      organism = 'euk',
                      run_mode = "starter")
            expect_is(check_khdel(sp, pattern = 'prosite'), 'ErResult')
            expect_message(check_khdel(sp, pattern = 'prosite'),
                           'Submitted sequences... 1')
            expect_message(
              check_khdel(sp, pattern = 'strict'),
              'Sequences with terminal ER retention signals detected... 0'
            )

            # check with fasta containing KDEL/HDEL motifs
            br <-
              readAAStringSet(system.file("extdata", "er_prot.fasta", 
                                          package = "SecretSanta"),
                              use.names = TRUE)
            inp <- setInfasta(inp, br)
            sp2 <- signalp(inp,
                           version = 2,
                           organism = 'euk',
                           run_mode = "starter")
            
            expect_message(
              check_khdel(sp2, pattern = 'prosite'),
              'Sequences with terminal ER retention signals detected... 2'
            )
            
            
          })