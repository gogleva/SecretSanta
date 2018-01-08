context("Piping tests")

test_that("workflows work",
          { # ---- This test takes a long time to run
            # ----- Minimal workflow: sp4 -> tmhmm -> check (K/H)DEL
            
            inp <- CBSResult()
            aa <- readAAStringSet(system.file("extdata", "small_prot.fasta", package = "SecretSanta"), use.names = TRUE)
            inp <- setInfasta(inp, aa)
            emp <- CBSResult()
            
            # ------- #Step1: signalp4
            s1_sp4 <- signalp(inp, version = 4, 'euk', run_mode = "starter")
            expect_is(s1_sp4, 'SignalpResult')
            
            # ------- #Step2: TMHMM with 0 or 1 TM domains allowed in mature peptide
            s2_tm0 <- tmhmm(s1_sp4, TM = 0)
            expect_is(s2_tm0, 'TMhmmResult')
            s2_tm1 <- tmhmm(s1_sp4, TM = 1)
            expect_is(s2_tm1, 'TMhmmResult')
            
            # ------- #Step3: check (K/H)DELs
            
            result0 <- check_khdel(s2_tm0, pattern = 'strict')
            result1 <- check_khdel(s2_tm1, pattern = 'strict')
            expect_is(result0, 'ErResult')
            expect_is(result1, 'ErResult')
            
            # ----- Stringent workflow: sp2 -> sp3 -> sp4 -> TMHMM -> wolf -> (K/H)DEL
            
            # ------- #Step1: signalp2
            s1_sp2 <- signalp(inp, version = 2, 'euk', run_mode = "starter",
                              legacy_method = 'hmm')
            expect_is(s1_sp2, 'SignalpResult')
            
            # ------- #Step2: signalp3
            
            s2_sp3 <- signalp(s1_sp2, version = 3, 'euk', run_mode = 'piper',
                              legacy_method = 'nn')
            expect_is(s2_sp3, 'SignalpResult')
            
            # ------- #Step3: signalp4
            
            s3_sp4 <- signalp(s2_sp3, version = 4, 'euk', run_mode = 'piper')
            expect_is(s3_sp4, 'SignalpResult')
            
            # ------- #Step4: TMHMM
            
            s4_tm <- tmhmm(s3_sp4, TM = 0)
            expect_is(s4_tm, 'TMhmmResult')
            
            # ------- #Step5: WoLFPsort
            
            s5_wo <- wolfpsort(s4_tm, organism = 'fungi', run_mode = 'piper' )
            expect_is(s5_wo, 'WolfResult')
            
            # ------- #Step6: Check C-terminal ER-retention signals:
            
            s6_result <- check_khdel(s5_wo, pattern = 'strict')
            expect_is(s6_result, 'CBSResult')
            
            
            # ----- Reverse workflow:  Wolf -> signalp4 -> TMHMM -> (K/H)DEL 
    
            # ------ #Step1: run WoLFPsort on the output:
            
            s1_wo <- wolfpsort(inp, organism = 'fungi', run_mode = 'starter')
            expect_is(s1_wo, 'WolfResult')
            
            # ------ #Step2: run signalp4:
            
            s2_sp4 <- signalp(s1_wo, organism = "euk", run_mode = 'piper', version = 4)
            
            # ------ #Step3: run TMHMM:
            
            s3_tm <- tmhmm(s2_sp4, TM = 0)
            
            # ------ #Step4: check (K/H)DEL
            
            s4_khdel <- check_khdel(s3_tm, pattern = 'strict')
            expect_is(s4_khdel, 'ErResult')
            
           
            # ----- Illegal workflows
            
            #------- #Start with TMHMM:
            
            expect_error(tmhmm(input_obj = inp, TM = 0),
                         'input_object does not belong to SignalpResult class')
            
            expect_error(tmhmm(input_obj = s4_khdel, TM = 0),
                         'input_object does not belong to SignalpResult class')
            
            expect_error(tmhmm(input_obj = s3_tm, paths = my_pa, TM = 1),
                         'input_object does not belong to SignalpResult class')
            
            # ----- Exhaustive input tests
            
            # Create missing input types:
            
            s1_sp3 <- signalp(inp, version = 3, 'euk', run_mode = "starter",
                              legacy_method = 'hmm')
            expect_is(s1_sp3, 'SignalpResult')
            
            s1_tp <- targetp(inp, network = 'N', run_mode = 'starter')
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
            # s2_wo - WolfPsort output, WolfResult object 
            
            # should be the same for all the tools, except TMHMM - this one can only pipe, so requires mature_fasta
            
            valid_inputs <- list(s1_sp2, s1_sp3, s1_sp4, s2_tm0, s1_wo, s1_tp)
            
            check_sp <- function(x, v, m) {
                result  <-  suppressMessages(
                                signalp(x, version = v, organism = 'euk',
                                        run_mode = m, legacy_method = 'hmm')
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
                                                    targetp(x, network = n,
                                                    run_mode = m)
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
            
            check_ER <- function(x, m) { 
                result  <-  suppressMessages(check_khdel(x, pattern = 'strict'))
                return(is(result, 'CBSResult')) 
            }
            
            #er_starters <- lapply(valid_inputs, check_ER)
            #expect_false(all(unlist(er_starters)))
            
            er_pipers <- lapply(valid_inputs, check_ER)
            expect_true(all(unlist(er_pipers)))
            
            
            tm_inputs <- list(s1_sp2, s1_sp3, s1_sp4)
            
            # ----- exhaustive unit tests for TMHMMm, should have mature fastas
            
            check_TMHMM <- function(x, t) {
                result <- suppressMessages(tmhmm(x, TM = t))
                return(is(result, 'CBSResult'))
            }
            
            tm_pipers <- lapply(tm_inputs, check_TMHMM, t = 0)
            expect_true(all(unlist(tm_pipers)))
            
            # ----  exhaustive unit tests for WolfPsort, should have out_fastas, so could be pipers only
            
            check_wolf <- function(x, o) { 
                result <- suppressMessages(wolfpsort(x, organism = o))
                return(is(result, 'CBSResult'))
            
            wolf_pipers <- lapply(valid_inputs, check_wolf, o = 'fungi')
            expect_true(all(unlist(wolf_pipers)))
            
            }
            })


