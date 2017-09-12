context("test helper functions for parallel sinalp")

test_that("combine_SignalpResult works as intended",
          {
            aa <- readAAStringSet(system.file("extdata",
                                              "sample_prot_100.fasta",
                                              package = "SecretSanta"))
            inp2 <- CBSResult(in_fasta = aa[1:10])
            inp3 <- CBSResult(in_fasta = aa[20:30])
            inp4 <- CBSResult(in_fasta = aa[40:50])
            
            sp1 <- signalp(
              input_obj = inp2,
              version = 2,
              organism = 'euk',
              run_mode = 'starter'
            )
            
            sp2 <- signalp(
              input_obj = inp3,
              version = 2,
              organism = 'euk',
              run_mode = 'starter'
            )
            sp3 <- suppressWarnings(signalp(
              input_obj = inp4,
              version = 2,
              organism = 'euk',
              run_mode = 'starter'
            ))
            
            obj <- list(sp1, sp2, sp3)
            combined_sp <- combine_SpResult(obj)
            
            expect_is(combined_sp, 'SignalpResult')
            
            expect_identical(sum(unlist(lapply(
              lapply(c(sp1, sp2, sp3), getInfasta), length
            ))),
            length(getInfasta(combined_sp)))
            
            expect_identical(sum(unlist(lapply(
              lapply(c(sp1, sp2, sp3), getOutfasta), length
            ))),
            length(getOutfasta(combined_sp)))
            
            expect_identical(sum(unlist(lapply(
              lapply(c(sp1, sp2, sp3), getSPtibble), nrow
            ))),
            nrow(getSPtibble(combined_sp)))
            
            expect_identical(sum(unlist(lapply(
              lapply(c(sp1, sp2, sp3), getMatfasta), length
            ))),
            length(getMatfasta(combined_sp)))
            
            alist <- c(inp2, inp3, inp4)
            expect_error(combine_SpResult(alist))
            
          })