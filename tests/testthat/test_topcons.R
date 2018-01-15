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
              
              # check object with empty out_fasta:
              
              expect_error(check_khdel(inp, pattern = 'strict'),
                           'the input object contains empty out_fasta slot'
              )
              
              # chek outputs of signalp
              sp <-
                  signalp(inp,
                          version = 2,
                          organism = 'euk',
                          run_mode = "starter",
                          legacy_method = 'hmm')
              
              
              
              expect_is(check_khdel(sp, pattern = 'prosite'), 'ErResult')
              expect_message(check_khdel(sp, pattern = 'prosite'),
                             'Submitted sequences... 1')
              expect_message(
                  check_khdel(sp, pattern = 'strict'),
                  'Sequences with terminal ER retention signals detected... 0'
              )
              
              # check with fasta containing KDEL/HDEL motifs
              br <-
                  readAAStringSet(system.file("extdata", "er_prot.fasta", 
                                              package = "SecretSanta"),
                                  use.names = TRUE)
              inp <- setInfasta(inp, br)
              sp2 <- signalp(inp,
                             version = 2,
                             organism = 'euk',
                             run_mode = "starter",
                             legacy_method = 'hmm')
              
              expect_message(
                  check_khdel(sp2, pattern = 'prosite'),
                  'Sequences with terminal ER retention signals detected... 2'
              )
              
              
          })