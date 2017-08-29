#' manage_paths function
#'
#' This function reads and stores pathways for external tools used in secertome
#' prediction pipeline and checks that all the external dependencies can be 
#' executed in principle, i.e produce correct help messages or process micro
#' fasta file. Required tools are available for academic users and could be
#' found here, for the instalation instructions see SecretSanta vignette.
#' \itemize{
#' \item WoLFPSORT - \url{https://github.com/fmaguire/WoLFPSort.git}
#' \item signalp2  - \url{http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp+2.0}
#' \item signalp3  - \url{http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp+3.0}
#' \item signalp4  - \url{http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp}
#' \item TMHMM     - \url{http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm}
#' \item targetP   - \url{http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?targetp}
#' }
#' @param check_installed if TRUE manage_paths will attempt to run external
#' dependencies, assuming that respective paths are attached to the $PATH
#' variable. Alteratively you can set check_installed to FALSE and supply 
#' path_file.
#' @param path_file  file paths to 2-column space-separated text file with 
#' listed paths for external dependencies if not sotred in $PATH;
#' \itemize{
#' \item first column should contain tool name;
#' \item second column should contain full path to the tool's executable;
#' \item for multiple versions of signalp use 'signalpV', 
#' where V is version number; 
#' }
#' @return tibble with verified file paths
#' @export
#' @examples
#' secret_paths <- manage_paths(system.file("extdata",
#'                                          "sample_paths",
#'                                          package="SecretSanta"))

manage_paths <- function(check_installed, path_file = NULL) {
  
  # here we assume that all the tools are acessible via $PATH:
  if (check_installed == TRUE) {
    message('checking dependencies acessible via $PATH')
    if (is.null(path_file)) {} else {
      message('path file provided, but not required')
    }
  }

  # alternatively, we will be dealing with paths supplied in a separate 
  # path_file, first check that all the supplied paths are valid
  if (check_installed == FALSE) {
    if (is.null(path_file)) {
      stop('Please supply path_file')
    } else {  
    # read path file in a tibble
    pp <- suppressMessages(readr::read_delim(path_file,
                                            delim = ' ',
                                            col_names = FALSE))
   
    # check that there are only 2 columns, i.e no extra spaces in tool names
    if (!(ncol(pp) == 2)) {
    stop('Please ensure that there are no spaces in the tool names.')}
   
    names(pp) <- c("tool", "path")
    # check that all supplied paths exist
    pp$status <- file.exists(pp$path)
   
    if (all(pp$status)) { 
    message('All paths are valid')
    # convert all the tool names to lower case to avoid confusion
    pp$tool <- tolower(pp$tool)
     
    } else {
    message('Supplied file path does not exist.')
    message(sapply(pp %>%
                  dplyr::filter_(~status == FALSE) %>%
                  dplyr::select_(~path), paste, '\n'))
    message('Please check that supplied paths are correct and try again')
    }
   }
  }
 
  
  # check that all the dependencies are executable in principle,
  # i.e we are able to run -h or
  # process a small sample fasta file
   
  my_tools <- c('signalp2',
               'signalp3',
               'signalp4',
               'targetp',
               'tmhmm',
               'wolfpsort')
  
  # micro fasta file to test with all tools:
  test_fasta <- system.file("extdata", "small_prot.fasta",
                            package = "SecretSanta")
  
  # helper function to extract tool paths when provided in path_file
  get_paths <- function(tool_name) {pp %>% 
                                    filter_(~tool == tool_name) %>%
                                    select_(~path)}
  
  make_call <- function(tool) {
    
    if (check_installed == TRUE) {
      tool <- tool
    } else if (check_installed == FALSE) {
      tool <- get_paths(tool)
    }
  }

  # now we will start making calls and evaluating their outputs:
  
  # signalp2:
  sp2_call <- system(paste(make_call('signalp2'), '-t euk', test_fasta),
                     intern = TRUE)
  if (grepl('SignalP predictions', sp2_call[1])) {
    message('siganlp2 run completed')
  } else {
    message('signalp2 test run failed; check if it is installed correctly')
  }
  
  #signalp3:
  
  sp3_call <- system(paste(make_call('signalp3'), '-t euk', test_fasta),
                     intern = TRUE)
  if (grepl('SignalP 3.0 predictions', sp3_call[1])) {
    message('siganlp3 run completed')
  } else {
    message('signalp3 test run failed; check if it is installed correctly')
  }
  
  #signalp4:
 
  sp4_call <- system(paste(make_call('signalp4'), '-t euk', test_fasta),
                     intern = TRUE)
  if (grepl('SignalP-4', sp4_call[1])) {
    message('siganlp4 run completed')
  } else {
    message('signalp4 test run failed; check if it is installed correctly')
  }
  
  #targetp:
  
  tp_call <- system(paste(make_call('targetp'), '-P', test_fasta),
                     intern = TRUE)
  if (grepl('targetp v1.1 prediction results', tp_call[2])) {
    message('targetp run completed')
  } else {
    message('targetp test run failed; check if it is installed correctly')
  }
  
  # tmhmm:
 
  tm_call <- system(paste(make_call('tmhmm'), '--short', test_fasta),
                    intern = TRUE)
  if (grepl('ALI_PLTG_1\tlen=94\tExpAA=22.44', tm_call[1])) {
    message('tmhmm run completed')
  } else {
    message('tmmhmm test run failed; check if it is installed correctly')
  }
  
  # wolfpsort:
  
  wolf_call <- system(paste(make_call('runWolfPsortSummary'), 'fungi', '>',
                            test_fasta),
                    intern = TRUE)
  
  # wolf_call <- system(paste('runWolfPsortSummary', 'fungi', '<', test_fasta),
  #                     intern = TRUE)
  if (grepl('ALI_PLTG_1\tlen=94\tExpAA=22.44', tm_call[1])) {
    message('tmhmm run completed')
  } else {
    message('tmmhmm test run failed; check if it is installed correctly')
  }
  
}    

  # sp2_call <- system(paste('signalp2', '-t euk', test_fasta), intern = TRUE)
  # sp3_call <- system(paste('signalp3', '-t euk', test_fasta), intern = TRUE)
  # sp4_call <- system(paste('signalp4', '-t euk', test_fasta), intern = TRUE)
  # tp_call <- system(paste('targetp', '-P', test_fasta), intern = TRUE)
  # tm_call <- system(paste('tmhmm', '--short', test_fasta), intern = TRUE)
  # wolf_call <- system(paste('runWolfPsortSummary', 'fungi', '<', test_fasta),
  #                     intern = TRUE)
  # 
  # } else if (check_installed == FALSE) {
  # sp2_call <- system(paste(get_paths('signalp2'), '-t euk', test_fasta),
  #                    intern = TRUE)
  # sp3_call <- system(paste('signalp3', '-t euk', test_fasta), intern = TRUE)
  # sp4_call <- system(paste('signalp4', '-t euk', test_fasta), intern = TRUE)
  # tp_call <- system(paste('targetp', '-P', test_fasta), intern = TRUE)
  # tm_call <- system(paste('tmhmm', '--short', test_fasta), intern = TRUE)
  # wolf_call <- system(paste('runWolfPsortSummary', 'fungi', '<', test_fasta),
  #                       intern = TRUE)  
  # }
  
  
  

  #helper function to produce status messages:
  check_signalp <- function(signalp_call, check_tool){
    if (check_tool == 'signalp2' || check_tool == 'signalp3') {
     if (signalp_call == 
         'Usage: signalp -t euk|gram+|gram- [options] seqfile'){
       message(paste(check_tool,
                     'test run completed'))
     } else {
       message(paste(check_tool,
                     'test run failed; check if it is installed correctly'))
     }
    } else if (check_tool == 'signalp4') {
      if (signalp_call == 
          '  Description: Predict signal peptide and cleavage site.'){
        message(paste(check_tool,
                      'test run completed'))
      } else {
        message(paste(check_tool,
                      'test run failed; check if it is installed correctly'))
      }
    }
   }
   
  check_signalp(call_signalp('signalp2'), 'signalp2')
  check_signalp(call_signalp('signalp3'), 'signalp3')
  check_signalp(call_signalp('signalp4'), 'signalp4')
}

  # 
  # # targetp 
  # 
  # if (suppressWarnings(system(paste(get_paths('targetp'), '-P',
  #                                   system.file("extdata",
  #                                               "small_prot.fasta",
  #                                               package = "SecretSanta")),
  #                             intern = TRUE)[9]) ==
  #     'ALI_PLTG_1             94   0.294  0.175  0.018  0.738   _    3') {
  #  message('targetp test run completed')
  # } else  {
  #   message('targetp test run failed; check if targetp is installed correctly')
  # }
  # 
  # # tmhmm
  # tm_call <- system(paste(get_paths('tmhmm'),
  #              system.file("extdata",
  #                          "small_prot.fasta", package="SecretSanta"),
  #              '--short'),
  #               intern = TRUE)
  # 
  # if (all(grepl('PredHel', tm_call))) {
  #  message('tmhmm test run completed')
  # } else {
  #   message('tmhmm test run failed; check if it is installed correctly')
  # }
  # 
  # # wolfpsort
  # 
  # wolf_call <- system(paste(get_paths('wolfpsort'),
  #                           'fungi', 
  #                           '<',
  #                           system.file("extdata",
  #                                       "small_prot.fasta",
  #                                       package="SecretSanta")),
  #                           intern = TRUE)
  # 
  # if (grepl('# k used for kNN is: 27', wolf_call[1])) {
  #   message('wolfpsort test run completed')
  # } else {
  #   message('wolfpsort test run failed; check if it is installed correctly')
  # }
  # 
  # return(pp)
#}