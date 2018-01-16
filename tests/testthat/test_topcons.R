context("Check TOPCONS integration")

test_that("TOPCONS outputs are integrated seamlessly",
          {
              # prep inputs:
              aa <- readAAStringSet(
                  system.file("extdata",
                              "sample_prot_100.fasta",
                              package = "SecretSanta"),
                  use.names = TRUE
              )
              p_dir <- system.file("extdata",
                                   "rst_SVw4hG.zip",
                                   package = "SecretSanta")
              
              inp <- CBSResult(in_fasta = aa)
              
              # try to input AAstringSet object:
              expect_error(topcons(input_obj = aa,
                             parse_dir = p_dir,
                             topcons_mode = "WEB-server",
                             TM = 0,
                             SP = FALSE),
                           "input_object does not belong to CBSResult class")
              
              # try to input CBSResult object with empty out_fasta slot
              
              expect_error(topcons(input_obj = inp,
                                   parse_dir = p_dir,
                                   topcons_mode = 'WEB-server',
                                   TM = 0),
                                   'out_fasta attribute is empty')
              # wrong file path
              expect_error(topcons(input_obj = inp,
                                   parse_dir = paste(p_dir, '000', sep = ''),
                                   topcons_mode = 'WEB-server',
                                   TM = 0),
              'Please provide valid path to the zipped TOPCONS output')
              
              # create some meaningful CBSResult object
              sp <-
                  signalp(inp,
                          version = 2,
                          organism = 'euk',
                          run_mode = "starter",
                          legacy_method = 'hmm')
              
              tpc <- topcons(input_obj = sp, parse_dir = p_dir, topcons_mode = "WEB-server", TM = 0, SP = TRUE)
              
              expect_error(topcons(input_obj = sp,
                      parse_dir = p_dir,
                      topcons_mode = "WEB-server",
                      TM = 0,
                      SP = 1),
                      "SP argument must be logical")
              
              expect_warning(topcons(input_obj = sp,
                             parse_dir = p_dir,
                             topcons_mode = "WEB-server",
                             TM = 2,
                             SP = FALSE),
                             "Recommended TM threshold values for secreted peptides is 0")
              
              expect_is(tpc, 'TopconsResult')
              expect_true(nrow(getTOPtibble(tpc)) == 9)
              
              
              # try with a larger file

          })


