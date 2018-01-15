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
              
              inp <- CBSResult(in_fasta = aa[1:10])
              
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
                          version = 4,
                          organism = 'euk',
                          run_mode = "starter",
                          legacy_method = 'hmm')
              
              tpc <- topcons(input_obj = sp,
                      parse_dir = p_dir,
                      topcons_mode = "WEB-server",
                      TM = 0,
                      SP = FALSE)
              
              expect_is(tpc, 'TopconsResult')
              expect_true(nrow(getTOPtibble(tpc)) == 1)

          })