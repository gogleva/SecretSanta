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
#' variable. Alteratively ifyou can set check_installed to FALSE and supply 
#' path_file.
#' @param path_file  file paths to external dependencies in a 2-column space-separated text file;
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
#'                                          

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
    stop('Please check that supplied paths are correct and try again')
    }
   }
  }
  
  # check that all the dependencies are executable in principle,
  # i.e we are able to process a small sample fasta file

  # micro fasta file to test with all tools:
  test_fasta <- system.file("extdata", "small_prot.fasta",
                            package = "SecretSanta")
  
  # helper function to extract tool paths when provided in path_file
  get_paths <- function(tool_name) {pp %>% 
                                    filter_(~tool == tool_name) %>%
                                    select_(~path)}
  
  # helper function to generate success message
  ok_message <- function(tool_name) {message(paste(tool_name,
                                                'run completed'))}
  # helper function to generate failure message
  failure_message <- function(tool_name) {message(
                                          paste(tool_name,
                                                'test run failed'))}
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
    ok_message('signalp2')
  } else {
    failure_message('signalp2')
  }

  #signalp3:
  sp3_call <- system(paste(make_call('signalp3'), '-t euk', test_fasta),
                     intern = TRUE)
  if (grepl('SignalP 3.0 predictions', sp3_call[1])) {
    ok_message('signalp3')
  } else {
    failure_message('signalp3')
  }

  #signalp4:
  sp4_call <- system(paste(make_call('signalp4'), '-t euk', test_fasta),
                     intern = TRUE)
  if (grepl('SignalP-4', sp4_call[1])) {
    ok_message('signalp4')
  } else {
    failure_message('signalp4')
  }

  #targetp:
  tp_call <- system(paste(make_call('targetp'), '-P', test_fasta),
                     intern = TRUE)
  if (grepl('targetp v1.1 prediction results', tp_call[2])) {
    ok_message('targetp')
  } else {
    failure_message('targetp')
  }

  # tmhmm:
  tm_call <- system(paste(make_call('tmhmm'), '--short', test_fasta),
                    intern = TRUE)
  if (grepl('ALI_PLTG_1\tlen=94\tExpAA=22.44', tm_call[1])) {
    ok_message('tmhmm')
  } else {
    failure_message('tmhmm')
  }

  # wolfpsort:
  wolf_call <- system(paste(make_call('wolfpsort'), 'fungi', '<',
                            test_fasta),
                      intern = TRUE)

  if (grepl('# k used for kNN is: 27', wolf_call[1])) {
    ok_message('wolfpsort')
  } else {
    failure_message('wolfpsort')
  }

}
