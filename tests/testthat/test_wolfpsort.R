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
          er <- check_khdel(tm, run_mode = 'piper')
          er_native <- check_khdel(inp, run_mode = 'starter')
          
          # run tests
          expect_is(wolfpsort(step1_sp2, 'fungi'), 'WolfResult')
          expect_is(wolfpsort(tm, 'fungi'), 'WolfResult')
          expect_is(wolfpsort(step2_sp3, 'fungi'), 'WolfResult')
          expect_is(wolfpsort(step1_sp2, 'plant'), 'WolfResult')
          expect_is(wolfpsort(step1_sp2, 'animal'), 'WolfResult')
          expect_is(wolfpsort(er, 'animal'), 'WolfResult')
          expect_is(wolfpsort(er_native, 'animal'), 'WolfResult')
          
          expect_error(wolfpsort(aa[1:10], 'fungi'),
                       'input_object does not belong to CBSResult superclass')
          expect_error(wolfpsort(inp, 'fungi'),
                       'the input object contains empty out_fasta slot')
          expect_error(wolfpsort(step1_sp2, 'plants'))
          
          })