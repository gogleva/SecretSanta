context("Check wolfpsort")

test_that("wolfpsort is ok and handles invalid inputs correctly",
          {
          # prep the inputs
          aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", 
                                            package = "SecretSanta"),
                                use.names = TRUE)
          inp <- CBSResult(in_fasta = aa[1:10])
          step1_sp2 <- signalp(inp,
                               version = 2,
                               organism = 'euk',
                               run_mode = "starter")
          step2_sp3 <- signalp(step1_sp2,
                               version = 3,
                               organism = 'euk',
                               run_mode = 'piper')
          #frozen tests
          tm <- tmhmm(step1_sp2, TM = 1)
          er <- check_khdel(tm)

          # run tests
          expect_is(wolfpsort(inp, 'fungi', run_mode = 'starter'), 'WolfResult')
          expect_is(wolfpsort(step1_sp2, 'fungi', run_mode = 'starter'), 'WolfResult')
          expect_is(wolfpsort(step1_sp2, 'fungi', run_mode = 'piper'), 'WolfResult')
          expect_is(wolfpsort(tm, 'fungi', run_mode = 'piper'), 'WolfResult')
          expect_is(wolfpsort(step2_sp3, 'fungi', run_mode = 'piper'), 'WolfResult')
          expect_is(wolfpsort(step1_sp2, 'plant', run_mode = 'piper'), 'WolfResult')
          expect_is(wolfpsort(step1_sp2, 'animal', run_mode = 'piper'), 'WolfResult')
          expect_is(wolfpsort(er, 'animal', run_mode = 'piper'), 'WolfResult')

          
          expect_error(wolfpsort(aa[1:10], 'fungi', run_mode = 'piper'),
                       'input_object does not belong to CBSResult superclass')
          
          expect_error(wolfpsort(step1_sp2, 'plants', run_mode = 'piper'))
          
          })