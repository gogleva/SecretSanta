context("test helper functions for parallel sinalp")

test_that("test split_XStringSet",
          {
            large_aa <- readAAStringSet(system.file("extdata", "Ppalm_prot_ALI_PLTG.fasta", package = "SecretSanta"))
            large_inp <- CBSResult(in_fasta = large_aa)
            res <- split_XStringSet(large_aa, 1000)
            expect_is(res, 'list')
            expect_error(split_XStringSet(large_aa, 50000), 'Chunk size exceeds total seq number')
            expect_error(split_XStringSet(large_inp, 100), 'Input string_set does not belong to XStringSet class')
          
          })

test_that("combine_SignalpResult combines",
          {
            inp2 <- CBSResult(
              in_fasta = readAAStringSet(system.file("extdata", "tail_prot.fasta", package = "SecretSanta")))
            inp3 <- CBSResult(
              in_fasta = readAAStringSet(system.file("extdata", "tail2_prot.fasta", package = "SecretSanta")))
            inp4 <- CBSResult(
              in_fasta = readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta")))

            sp1 <- signalp(input_obj = inp2, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa)
            sp2 <- signalp(input_obj = inp3, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa)
            sp3 <- signalp(input_obj = inp4, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa)
            obj <- list(sp1, sp2, sp3)
            combined_sp <- combine_SignalpResult(obj)
            
            expect_identical(sum(unlist(lapply(lapply(c(sp1, sp2, sp3), getInfasta), length))), 
                             length(getInfasta(combined_sp)))
            expect_identical(sum(unlist(lapply(lapply(c(sp1, sp2, sp3), getOutfasta), length))), 
                             length(getOutfasta(combined_sp)))
            expect_identical(sum(unlist(lapply(lapply(c(sp1, sp2, sp3), getSPtibble), nrow))), 
                             length(getSPtibble(combined_sp)))
            expect_identical(sum(unlist(lapply(lapply(c(sp1, sp2, sp3), getMatfasta), length))), 
                             length(getMatfasta(combined_sp)))
                                                 
          })

