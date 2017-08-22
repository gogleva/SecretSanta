context("Check tmhmm")

test_that("tmhmm outputs are correctly parsed for system calls",
          {
          # prep inut object:
          
          my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
          inp <- SignalpResult()
          aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"),
                                use.names = TRUE)
          inp <- setInfasta(inp, aa)
          s1_sp2 <- signalp(inp, version = 2, 'euk', run_mode = "starter", paths = my_pa)
            
          # run tmhmm on the output of signalp step  
          expect_is(tmhmm(s1_sp2, paths = my_pa, TM = 1), 'TMhmmResult')
          
          # run tmhmm on XStringSetObject:
          expect_error(tmhmm(aa, paths = my_pa, TM = 0), 'input_object does not belong to SignalpResult class')
          
          # run tmhmm on CBSresult object with empty mature_fasta slot:
          expect_error(tmhmm(inp, paths = my_pa, TM = 0), 'the input object contains an empty mature_fasta slot')
          
          # run tmhmm on object without mature_fasta slot
          cb <- CBSResult()
          cb <- setInfasta(cb, aa)
          expect_error(tmhmm(cb, paths = my_pa, TM = 0), 'input_object does not belong to SignalpResult class')
          
          })