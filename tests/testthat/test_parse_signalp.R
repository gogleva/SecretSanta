context("Check signalp parsing")

test_that("signalp outputs are correctly parsed for system calls",
          {
            s_fasta <- system.file("extdata", "sample_prot.fasta", package = "SecretSanta") 
            secret_paths <- manage_paths(system.file("extdata", "sample_paths", package="SecretSanta"))
            sp2_path <- secret_paths %>% filter(tool == 'signalp2') %>% select(path)
            con <- system(paste(sp2_path, '-t euk', s_fasta), intern = TRUE)
            parse_sp_system <- parse_signalp(input = con, input_type = "system_call")
            expect_is(parse_sp_system, 'tbl')
            expect_equal(ncol(parse_sp_system), 9)
            expect_true(all(parse_sp_system$Prediction == 'Signal peptide'))
            
            #check that all the fields are present:
            fields <- c("gene_id", "Cmax", "Cpos", "Ymax", "Ypos",
                         "Smax", "Spos", "Smean", "Prediction")
            expect_equal(names(parse_sp_system), fields)
            
            #check that gene_ids do not contain extra spaces:
            expect_false(' ' %in% parse_sp_system$gene_id)
            })

test_that("signalp outputs are correctly parsed for files",
          {
            s_path <- system.file("extdata", "sample_prot_signalp2_out", package = "SecretSanta") 
            parse_sp_path <- parse_signalp(input = s_path, input_type = "path")
            expect_is(parse_sp_path, 'tbl')
            expect_equal(ncol(parse_sp_path), 9)
            expect_true(all(parse_sp_path$Prediction == 'Signal peptide'))
            
            #check that all the fields are present
            fields <- c("gene_id", "Cmax", "Cpos", "Ymax", "Ypos",
                        "Smax", "Spos", "Smean", "Prediction")
            expect_equal(names(parse_sp_path), fields)
            
            #check that gene_ids do not contain extra spaces:
            expect_false(' ' %in% parse_sp_path$gene_id)
            })