test_that("test split_XStringSet",
          { aa_10K <- readAAStringSet("~/anna/Labjournal/SecretSanta_external/test_fastas/large_10K.fasta")
            large_inp <- CBSResult(in_fasta = aa_10K)
            res <- split_XStringSet(aa_10K, 1000)
            expect_is(res, 'list')
            expect_error(split_XStringSet(aa_10K, 50000), 'Chunk size exceeds total seq number')
            expect_error(split_XStringSet(large_inp, 100), 'Input string_set does not belong to XStringSet class')
            
          })

test_that("signalp processes large inputs correctly",
          {          
            aa_10K <- readAAStringSet("~/anna/Labjournal/SecretSanta_external/test_fastas/large_10K.fasta")
            sp2_10K <- signalp(CBSResult(in_fasta = aa_10K[1:1000]),
                               organism = 'euk',
                               version = 2,
                               run_mode = 'starter',
                               truncate = T)
            
            rescued <- m_slicer(sp2_10K,  # signalp2 output
                                run_mode = 'rescue',
                                min_len = 100)
            
            sp2_rinp <- signalp(CBSResult(in_fasta = rescued),
                               organism = 'euk',
                               version = 2,
                               run_mode = 'starter')
            expect_is(sp2_rinp, "SignalpResult")
          }
          )
