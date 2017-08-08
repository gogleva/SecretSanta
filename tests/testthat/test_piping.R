context("Piping tests")

test_that("workflows work",
          {
            # ----- Minimal workflow: sp4 -> tmhmm -> check (K/H)DEL
            
            my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
            inp <- CBSResult()
            aa <- readAAStringSet(system.file("extdata", "sample_prot.fasta", package = "SecretSanta"),
                                  use.names = TRUE)
            inp <- setInfasta(inp, aa)
            emp <- CBSResult()
            
            # ------- #Step1: signalp4
            s1_sp4 <- signalp(inp, version = 4, 'euk', run_mode = "starter", paths = my_pa)
            expect_is(s1_sp4, 'SignalpResult')
            
            # ------- #Step2: TMHMM with 0 or 1 TM domains allowed in mature peptide
            s2_tm0 <- tmhmm(s1_sp4, paths = my_pa, TM = 0)
            expect_is(s2_tm0, 'TMhmmResult')
            s2_tm1 <- tmhmm(s1_sp4, paths = my_pa, TM = 1)
            expect_is(s2_tm1, 'TMhmmResult')
            
            # ------- #Step3: check (K/H)DELs
            
            result0 <- check_khdel(s2_tm0, run_mode = 'piper')
            result1 <- check_khdel(s2_tm1, run_mode = 'piper')
            expect_is(result0, 'ErResult')
            expect_is(result1, 'ErResult')
            
            # ----- Stringent workflow: sp2 -> sp3 -> sp4 -> TMHMM -> wolf -> (K/H)DEL
            
            # ------- #Step1: signalp2
            s1_sp2 <- signalp(inp, version = 2, 'euk', run_mode = "starter", paths = my_pa)
            expect_is(s1_sp2, 'SignalpResult')
            
            # ------- #Step2: signalp3
            
            s2_sp3 <- signalp(s1_sp2, version = 3, 'euk', run_mode = 'piper', paths = my_pa)
            expect_is(s2_sp3, 'SignalpResult')
            
            # ------- #Step3: signalp4
            
            s3_sp4 <- signalp(s2_sp3, version = 4, 'euk', run_mode = 'piper', paths = my_pa)
            expect_is(s3_sp4, 'SignalpResult')
            
            # ------- #Step4: TMHMM
            
            s4_tm <- tmhmm(s3_sp4, paths = my_pa, TM = 0)
            expect_is(s4_tm, 'TMhmmResult')
            
            # ------- #Step5: WoLFPsort
            
            s5_wo <- wolfpsort(s4_tm, organism = 'fungi', paths = my_pa)
            expect_is(s5_wo, 'WolfResult')
            
            # ------- #Step6: Check C-terminal ER-retention signals:
            
            s6_result <- check_khdel(s5_wo, run_mode = 'piper')
            expect_is(s6_result, 'CBSResult')
            
            
            # ----- Reverse workflow: (K/H)DEL -> Wolf -> signalp4 -> TMHMM
            
            # ------ #Step1: remove sequences with N-terminal ER-retention signals:
            
            s1_er <- check_khdel(inp, run_mode = 'starter')
            expect_is(s1_er, c('CBSResult', 'ErResult'))
            
            # ------ #Step2: run WoLFPsort on the output:
            
            s2_wo <- wolfpsort(s1_er, organism = 'fungi', paths = my_pa)
            expect_is(s2_wo, 'WolfResult')
            
            
            # ----- Illegal workflows
            
            #------- #Start with TMHMM:
            
            expect_error(tmhmm(input_obj = inp, paths = my_pa, TM = 0),
                         'the input object does not contain mature_fasta slot')
            
            expect_error(tmhmm(input_obj = s1_er, paths = my_pa, TM = 0),
                         'the input object does not contain mature_fasta slot')
            
            expect_error(tmhmm(input_obj = s2_wo, paths = my_pa, TM = 1),
                         'the input object does not contain mature_fasta slot')
            
            # ----- Exhaustive input tests
            
            # Create missing input types:
            
            s1_sp3 <- signalp(inp, version = 3, 'euk', run_mode = "starter", paths = my_pa)
            expect_is(s1_sp3, 'SignalpResult')
            
            s1_tp <- targetp(inp, network_type = 'N', run_mode = 'starter', paths = my_pa)
            expect_is(s1_tp, 'TargetpResult')
            
            # Inputs:
            
            # starter inputs, tested in the respective tool unit tests
            # aa - AAStringSet object - tested in signalp unit test
            # inp - CBSResult object with filled input_fasta - tested in signlap unit test
            # emp - empty CBSResult object - tested in signalp unit test
            
            # Outputs from other tools:
            # s1_sp2 - signalp2 output, SignalpResult object
            # s1_sp3 - signalp3 output, SignalpResult object
            # s1_sp4 - signalp4 output, SignalpResult object
            # s1_tp - targetp output, SignalpResult object
            # s2_tm0 - TMHMM output, TMhmmResult object # === PIPER
            # s1_er - check K/HDEL output, ErResult object
            # s2_wo - WolfPsort output, WolfResult object # === PIPER
            
            # should be the same for all the tools, except TMHMM - this one can only pipe, so requires mature_fasta
            
            valid_inputs <- list(s1_sp2, s1_sp3, s1_sp4, s2_tm0, s1_er, s2_wo, s1_tp)
            
            check_sp <- function(x, v, m) { result  <-  suppressMessages(
                                                        signalp(x, version = v, organism_type = 'euk',
                                                        run_mode = m, paths = my_pa)
                                                        )
                                            return(is(result, 'CBSResult')) 
                                     } 
            # ----- Exhaustive input tests for signalp2
            
            # starters
            sp2_starters <- lapply(valid_inputs, check_sp, v = 2, m = 'starter')
            expect_true(all(unlist(sp2_starters)))
            
            # pipers
            sp2_pipers <- lapply(valid_inputs, check_sp, v = 2, m = 'piper')
            expect_true(all(unlist(sp2_pipers)))
            
   
            # ----- Exhaustive input tests for signalp3
            
            # starters
            sp3_starters <- lapply(valid_inputs, check_sp, v = 3, m = 'starter')
            expect_true(all(unlist(sp3_starters)))
            
            # pipers
            sp3_pipers <- lapply(valid_inputs, check_sp, v = 3, m = 'piper')
            expect_true(all(unlist(sp3_pipers)))
            
            
            # ----- exhaustive input tests for signalp4
            
            # starters
            sp4_starters <- lapply(valid_inputs, check_sp, v = 4, m = 'starter')
            expect_true(all(unlist(sp4_starters)))
            
            # pipers
            sp4_pipers <- lapply(valid_inputs, check_sp, v = 4, m = 'piper')
            expect_true(all(unlist(sp4_pipers)))
            
            
            # ----- exhaustive input tests for targetp
            
            check_tp <- function(x, n, m) { result  <-  suppressMessages(
                                                    targetp(x, network_type = n,
                                                    run_mode = m, paths = my_pa)
                                                    )
                                            return(is(result, 'CBSResult')) 
                                           } 
            # starters:
            tp_starters_N <- lapply(valid_inputs, check_tp, n = 'N', m = 'starter')
            expect_true(all(unlist(tp_starters_N)))
            
            tp_starters_P <- lapply(valid_inputs, check_tp, n = 'P', m = 'starter')
            expect_true(all(unlist(tp_starters_P)))
            
            # pipers:
            
            tp_pipers_N <- lapply(valid_inputs, check_tp, n = 'N', m = 'piper')
            expect_true(all(unlist(tp_pipers_N)))
            
            tp_pipers_P <- lapply(valid_inputs, check_tp, n = 'P', m = 'piper')
            expect_true(all(unlist(tp_pipers_P)))
            
            
            
            # ----- exhaustive input tests for (K/H)DEL check:
            
            check_ER <- function(x, m) { result  <-  suppressMessages(check_khdel(x, run_mode = m))
                                      return(is(result, 'CBSResult')) 
            }
            
            er_starters <- lapply(valid_inputs, check_ER, m = 'starter')
            expect_true(all(unlist(er_starters)))
            
            er_pipers <- lapply(valid_inputs, check_ER, m = 'piper')
            expect_true(all(unlist(er_pipers)))
            
            })


