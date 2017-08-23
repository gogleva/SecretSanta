context("Check ER-retention motif detection")

test_that("terminal KHDEL/HDEL motifs are detected",
          {
            # prep inputs:
            inp <- CBSResult()
            aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), 
                                  use.names = TRUE)
            inp <- setInfasta(inp, aa[1:10])
            
            # check valid inputs, starter mode on a CBSResult object:
            expect_is(check_khdel(inp, run_mode = 'starter'), 'ErResult')
            
            
            # check object with empty out_fasta in piper run_mode:
            expect_error(check_khdel(inp, run_mode = 'piper'), 
                        'query fasta is empty, please ensure you are using correct run_mode')
            
            # chek outputs of signalp
            my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
            sp <- signalp(inp, version = 2, 'euk', run_mode = "starter", paths = my_pa)
            expect_is(check_khdel(sp, run_mode = 'piper'), 'ErResult')
            expect_message(check_khdel(sp, run_mode = 'piper'),
                           'Submitted sequences... 8')
            expect_message(check_khdel(sp, run_mode = 'piper'),
                           'Sequences with terminal ER retention signals detected... 0')
            
            # check outputs of TMHMM
            tm <- tmhmm(sp, paths = my_pa, TM = 1)
            expect_is(check_khdel(tm, run_mode = 'piper'), 'ErResult')
            
            # check with fasta containing KDEL/HDEL motifs
            br <- readAAStringSet(system.file("extdata", "er_prot.fasta", package = "SecretSanta"), 
                                  use.names = TRUE)
            inp <- setInfasta(inp, br)
            expect_is(check_khdel(inp, run_mode = 'starter'),  'ErResult')
            expect_message(check_khdel(inp, run_mode = 'starter'),
                           'Submitted sequences... 5')
            expect_message(check_khdel(inp, run_mode = 'starter'),
                           'Sequences with terminal ER retention signals detected... 2')
            

          })