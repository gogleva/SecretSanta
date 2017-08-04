context("Piping tests")

test_that("Minimal workflow works",
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
            
        })


test_that("Stringent workflow works",
          {
            
          })
