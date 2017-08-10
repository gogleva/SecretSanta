context("test helper functions for parallel sinalp")

test_that("test split_XStringSet",
          {
            large_aa <- readAAStringSet(system.file("extdata", "Ppalm_prot_ALI_PLTG.fasta", package = "SecretSanta"))
            large_inp <- CBSResult(in_fasta = large_aa)
            res <- split_XStringSet(large_aa, 1000)
            expect_is(res, 'list')
            expect_error(split_XStringSet(large_aa, 50000), 'Chunk size exceeds total seq number')
            expect_error(split_XStringSet(large_inp, 100))
          
          })