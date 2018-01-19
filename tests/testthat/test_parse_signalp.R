context("Check signalp parsing")

test_that("signalp outputs are correctly parsed for system calls",
          {
            s_fasta <-
              system.file("extdata", "sample_prot_100.fasta",
                          package = "SecretSanta")
            con_hmm_v2 <- system(paste('signalp2 -t euk -f short -m hmm -trunc 70', s_fasta), intern = TRUE)
            con_nn_v2 <- system(paste('signalp2 -t euk -f short -m nn -trunc 70', s_fasta), intern = TRUE)
            
            con_hmm_v3 <- system(paste('signalp3 -t euk -f short -m hmm -trunc 70', s_fasta), intern = TRUE)
            con_nn_v3 <- system(paste('signalp3 -t euk -f short -m nn -trunc 70', s_fasta), intern = TRUE)
            
            path_hmm_2 <- system.file("extdata",
                                      "sample_prot_sp2_hmm_out",
                                      package = "SecretSanta")
            
            path_nn_2 <- system.file("extdata",
                                     "sample_prot_sp2_nn_out",
                                     package = "SecretSanta")    
            
            path_hmm_3 <- system.file("extdata",
                                      "sample_prot_sp3_hmm_out",
                                      package = "SecretSanta")
            
            path_nn_3 <- system.file("extdata",
                                     "sample_prot_sp3_nn_out",
                                     package = "SecretSanta")
            
            test_v_con <- function(version, hmm_in, nn_in, i_type) {
            
            parse_hmm_system <-
              parse_signalp(input = hmm_in, input_type = i_type,
                            method = 'hmm', version = version)
            expect_is(parse_hmm_system, 'tbl')
            expect_equal(ncol(parse_hmm_system), 9)
            expect_true(all(parse_hmm_system$Prediction == 'Signal peptide'))

            expect_error(
                parse_signalp(input = nn_in, input_type = i_type,
                              method = 'nn'),
                'Missing argument: version.')
            
            if (version == 2) {
            expect_error(
                parse_signalp(input = nn_in, input_type = i_type,
                              method = 'nn', version = version),
                'Please provide source_fasta.')
            }
            
            parse_nn_system <-
                parse_signalp(input = nn_in, input_type = i_type,
                              method = 'nn', version = version, source_fasta = s_fasta)
            expect_is(parse_nn_system, 'tbl')
            expect_equal(ncol(parse_nn_system), 9)
            expect_true(all(parse_nn_system$Prediction == 'Signal peptide'))
            
            #check that all the fields are present:
            fields <- c("gene_id",
                        "Cmax",
                        "Cpos",
                        "Ymax",
                        "Ypos",
                        "Smax",
                        "Spos",
                        "Smean",
                        "Prediction")
            
            expect_equal(names(parse_hmm_system), fields)
            expect_equal(names(parse_nn_system), fields)
            
            #check that gene_ids do not contain extra spaces:
            expect_false(' ' %in% parse_hmm_system$gene_id)
            expect_false(' ' %in% parse_nn_system$gene_id)
            }
            
            # test system calls:
            # test signalp2
            test_v_con(2, con_hmm_v2, con_nn_v2, 'system_call')
            
            # test signalp3
            test_v_con(3, con_hmm_v3, con_nn_v3, 'system_call')
            
            # test parsing from file:
            # test signalp2
            test_v_con(2, path_hmm_2, path_nn_2, 'path')
            
            # test signalp3
            test_v_con(3, path_hmm_3, path_nn_3, 'path')
            
         })