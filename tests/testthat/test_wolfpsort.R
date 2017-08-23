context("Check wolfpsort")

test_that("wolfpsort outputs correct objects and handles invalid inputs correctly",
          {
          # prep the inputs
            
          my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
          inp <- SignalpResult()
          aa <- readAAStringSet(system.file("extdata", "sample_prot.fasta", package = "SecretSanta"), use.names = TRUE)
          inp <- setInfasta(inp, aa[1:10])
          step1_sp2 <- signalp(inp, version = 2, 'euk', run_mode = "starter", paths = my_pa)
          step2_sp3 <- signalp(step1_sp2, version = 3, 'euk', run_mode = 'piper', paths = my_pa)
          tm <- tmhmm(step1_sp2, paths = my_pa, TM = 1)
          er <- check_khdel(tm, run_mode = 'piper')
          er_native <- check_khdel(inp, run_mode = 'starter')
          
          # run tests
          expect_is(wolfpsort(step1_sp2, 'fungi', my_pa), 'WolfResult')
          expect_is(wolfpsort(tm, 'fungi', my_pa), 'WolfResult')
          expect_is(wolfpsort(step2_sp3, 'fungi', my_pa), 'WolfResult')
          expect_is(wolfpsort(step1_sp2, 'plant', my_pa), 'WolfResult')
          expect_is(wolfpsort(step1_sp2, 'animal', my_pa), 'WolfResult')
          expect_is(wolfpsort(er, 'animal', my_pa), 'WolfResult')
          expect_is(wolfpsort(er_native, 'animal', my_pa), 'WolfResult')
          
          expect_error(wolfpsort(aa, 'fungi', my_pa),
                       'input_object does not belong to CBSResult superclass')
          expect_error(wolfpsort(inp, 'fungi', my_pa),
                       'the input object contains empty out_fasta slot')
          expect_error(wolfpsort(step1_sp2, 'plants', my_pa),
                       'input organism is not allowed or does not exist')
          
          })