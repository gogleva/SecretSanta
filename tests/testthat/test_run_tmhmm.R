context("Check tmhmm")

test_that("tmhmm outputs are correctly parsed for system calls",
          {
          # prep inut object:
          aa <- readAAStringSet(system.file("extdata",
                                            "sample_prot_100.fasta",
                                            package = "SecretSanta"))
          inp <- CBSResult(in_fasta = aa[1:10])
          s1_sp2 <- signalp(inp,
                            version = 2, 
                            organism = 'euk',
                            run_mode = "starter",
                            legacy_method = 'hmm')
            
          # run tmhmm on the output of signalp step  
          expect_is(tmhmm(s1_sp2, TM = 1), 'TMhmmResult')
          
          # run tmhmm on XStringSetObject:
          expect_error(tmhmm(aa, TM = 0),
                       'Input object does not belong to SignalpResult class.')
          
          # run tmhmm on CBSresult object with empty mature_fasta slot:
          expect_error(tmhmm(inp, TM = 0),
                       'Input object does not belong to SignalpResult class.')
          
          # run tmhmm on object without mature_fasta slot
          cb <- CBSResult()
          cb <- setInfasta(cb, aa)
          expect_error(tmhmm(cb, TM = 0),
                       'Input object does not belong to SignalpResult class.')
          
          })