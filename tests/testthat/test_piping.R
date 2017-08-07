context("Piping tests")

test_that("workflows work",
          {
            # ----- Minimal workflow: sp4 -> tmhmm -> check (K/H)DEL
            
            my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
            inp <- SignalpResult()
            aa <- readAAStringSet(system.file("extdata", "sample_prot.fasta", package = "SecretSanta"),
                                  use.names = TRUE)
            inp <- setInfasta(inp, aa)
            
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
            
     #       expect_error(tmhmm(input_obj = inp, paths = my_pa, TM = 0),
      #                   'the input object contains an empty mature_fasta slot') => wrong error
            expect_error(tmhmm(input_obj = s1_er, paths = my_pa, TM = 0),
                         'the input object does not contain mature_fasta slot')
            expect_error(tmhmm(input_obj = s2_wo, paths = my_pa, TM = 1),
                         'the input object does not contain mature_fasta slot')
      
            })


