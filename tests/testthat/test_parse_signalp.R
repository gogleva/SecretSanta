context("Check signalp parsing")

test_that("signalp outputs are correctly parsed for system calls",
          {
            s_fasta <-
              system.file("extdata", "sample_prot_100.fasta",
                          package = "SecretSanta")
            con_hmm_v2 <- system(paste('signalp2 -t euk -f short -m hmm -trunc 70', s_fasta), intern = TRUE)
            con_nn_v2 <- system(paste('signalp2 -t euk -f short -m nn -trunc 70', s_fasta), intern = TRUE)
            
            parse_hmm2_system <-
              parse_signalp(input = con_hmm_v2, input_type = "system_call",
                            method = 'hmm', version = 2)
            expect_is(parse_hmm2_system, 'tbl')
            expect_equal(ncol(parse_hmm2_system), 10)
            expect_true(all(parse_hmm2_system$Prediction == 'Signal peptide'))

            expect_error(
                parse_signalp(input = con_nn_v2, input_type = "system_call",
                              method = 'nn'),
                'missing argument: version')
            
            expect_error(
                parse_signalp(input = con_nn_v2, input_type = "system_call",
                              method = 'nn', version = 2),
                'please provide source_fasta')
            
            parse_nn2_system <-
                parse_signalp(input = con_nn_v2, input_type = "system_call",
                              method = 'nn', version = 2, source_fasta = s_fasta)
            ### to do: continue with tests
            
            expect_is(parse_nn2_system, 'tbl')
            expect_equal(ncol(parse_nn2_system), 10)
            expect_true(all(parse_nn2_system$Prediction == 'Signal peptide'))
            
            #check that all the fields are present:
            fields <- c("gene_id",
                        "Cmax",
                        "Cpos",
                        "Ymax",
                        "Ypos",
                        "Smax",
                        "Spos",
                        "Smean",
                        "Prediction",
                        "Prediction_YN")
            expect_equal(names(parse_hmm2_system), fields)
            
            #check that gene_ids do not contain extra spaces:
            expect_false(' ' %in% parse_hmm2_system$gene_id)
          })

test_that("signalp outputs are correctly parsed for files",
          {
            s_path <-
              system.file("extdata", "sample_prot_signalp2_out2",
                          package = "SecretSanta")
            parse_sp_path <-
              parse_signalp(input = s_path, input_type = "path")
            expect_is(parse_sp_path, 'tbl')
            expect_equal(ncol(parse_sp_path), 9)
            expect_true(all(parse_sp_path$Prediction == 'Signal peptide'))
            
            #check that all the fields are present
            fields <- c("gene_id",
                        "Cmax",
                        "Cpos",
                        "Ymax",
                        "Ypos",
                        "Smax",
                        "Spos",
                        "Smean",
                        "Prediction")
            expect_equal(names(parse_sp_path), fields)
            
            #check that gene_ids do not contain extra spaces:
            expect_false(' ' %in% parse_sp_path$gene_id)
          })