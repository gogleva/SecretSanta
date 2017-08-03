context("Check ER-retention motif detection")

test_that("terminal KHDEL/HDEL motifs are detected",
          {
            # prep inputs:
            inp <- SignalpResult()
            aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), 
                                  use.names = TRUE)
            inp <- setInfasta(inp, aa)
            
            # check valid inputs, starter mode on a CBSResult object:
            expect_is(check_khdel(inp, run_mode = 'starter'), 'ErResult')
            
            # check object with empty out_fasta in piper run_mode:
            expect_error(check_khdel(inp, run_mode = 'piper'), 
                        'query fasta is empty, please ensure you are using correct run_mode')
            
            # check ER retention signal in the signalp output, 'starter' mode
            #et_sp <- check_khdel(step1_sp2, run_mode = 'starter')
            
            # check ER retention signal in the signalp output, 'piper' mode
            #et_piper <- check_khdel(step1_sp2, run_mode = 'piper')  
          })