context("Check ask_uniprot")

test_that("ask_uniprot works",
          {
           
           # OK ids:         
           id_list <- c('P39864', 'D0N4E2', 'Q5BUB4', 'D0N381', 'B1NNT7', 'D0NP26')     
           result_ok <- ask_uniprot(id_list)
           expect_is(result_ok, 'tbl_df')
           expect_equal(id_list, result_ok$UniprotID)
           expect_equal(dim(result_ok), c(6, 3))
           
           # some bad/non-existent ids:
           bad_list <- c('P39864', 'D0N4E2', 'Q5BUB4', 'D0N381', 'B1NNT7', 'D0NP26999')
           expect_warning(ask_uniprot(bad_list))
           result_404 <- ask_uniprot(bad_list)
           expect_is(result_404, 'tbl_df')
           expect_equal(bad_list, result_404$UniprotID)
           expect_equal(dim(result_404), c(6, 3))
           expect_identical(result_404$Subcellular.Location[6], 
                            "Client error: (404) Not Found")
           })


