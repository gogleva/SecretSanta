context("Check signalp parsing")

test_that("signalp outputs are correctly parsed for system calls",
          {
            s_file <- system.file("extdata", "sample_paths", package = "SecretSanta")
            con <- system(paste(sp2_path, '-t euk', s_fasta), intern = TRUE)
            parse_sp_system <- parse_signalp(input = con, input_type = "system_call")
            expect_is(parse_sp_system, 'tbl')
            expect_equal(ncol(parse_sp_system), 11)
            })

test_that("signalp outputs are correctly parsed for files",
          {
            s_path <- system.file("extdata", "sample_prot_signalp2_out", package = "SecretSanta") 
            parse_sp_path <- parse_signalp(input = s_path, input_type = "path")
            expect_is(parse_sp_path, 'tbl')
            expect_equal(ncol(parse_sp_path), 11)
          })