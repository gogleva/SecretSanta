context("test helper functions for parallel sinalp")

test_that("combine_SignalpResult works as intended",
          {
            inp2 <- CBSResult(
              in_fasta = readAAStringSet(system.file("extdata", "tail_prot.fasta", package = "SecretSanta")))
            inp3 <- CBSResult(
              in_fasta = readAAStringSet(system.file("extdata", "tail2_prot.fasta", package = "SecretSanta")))
            inp4 <- CBSResult(
              in_fasta = readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta")))
            
            my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
            sp1 <- signalp(input_obj = inp2, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa)
            sp2 <- suppressWarnings(
                   signalp(input_obj = inp3, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa)
                                   )
            sp3 <- signalp(input_obj = inp4, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa)
                   
            obj <- list(sp1, sp2, sp3)
            combined_sp <- combine_SignalpResult(obj)
            
            expect_is(combined_sp, 'SignalpResult')
            
            expect_identical(sum(unlist(lapply(lapply(c(sp1, sp2, sp3), getInfasta), length))), 
                             length(getInfasta(combined_sp)))
            expect_identical(sum(unlist(lapply(lapply(c(sp1, sp2, sp3), getOutfasta), length))), 
                             length(getOutfasta(combined_sp)))
            expect_identical(sum(unlist(lapply(lapply(c(sp1, sp2, sp3), getSPtibble), nrow))), 
                             length(getSPtibble(combined_sp)))
            expect_identical(sum(unlist(lapply(lapply(c(sp1, sp2, sp3), getMatfasta), length))), 
                             length(getMatfasta(combined_sp)))
            
            alist <- c(inp2, inp3, inp4)
            expect_error(combine_SignalpResult(alist), 
                         'Some objects from arguments list do not belong to SignalpResult class.' )
            
            # test combine CBS
            
            combined_cbs <- combine_CBSResult(inp2, inp3, inp4)
            expect_is(combined_cbs, 'CBSResult')
            
            expect_identical(sum(unlist(lapply(lapply(c(inp2, inp3, inp4), getInfasta), length))), 
                             length(getInfasta(combined_cbs)))
            expect_identical(sum(unlist(lapply(lapply(c(inp2, inp3, inp4), getOutfasta), length))), 
                             length(getOutfasta(combined_cbs)))
            
            expect_warning(combine_CBSResult(sp1, sp2, sp3), 'Only in_fasta and out_fasta slots will be combined')
                                                 
          })