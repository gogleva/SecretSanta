#' manage_paths function
#'
#' This function reads and stores pathways for external tools used in secertome prediction pipeline. 
#' Required tools could be found here:
#' \itemize{
#' \item WoLFPSORT - \url{https://github.com/fmaguire/WoLFPSort.git}
#' \item signalp2  - \url{http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp+2.0}
#' \item signalp3  - \url{http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp+3.0}
#' \item signalp4  - \url{http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp}
#' \item TMHMM     - \url{http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm}
#' \item targetP   - \url{http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?targetp}
#' }
#' @param path_file  2-column space-separated file with listed paths for external dependencies;
#' \itemize{
#' \item first column should contain tool name;
#' \item second column should contain full path to the tool's executable;
#' \item for multiple versions of signalp use 'signalpV', where V is version number; 
#' }
#' @return tibble with verified file paths
#' @export
#' @examples
#' secret_paths <- manage_paths(system.file("extdata", "sample_paths", package="SecretSanta"))

manage_paths <- function(path_file) {
  # read path file in a tibble
  pp <- suppressMessages(readr::read_delim(path_file, delim = ' ', col_names = FALSE))
  
  # check that there are only 2 columns, i.e no extra spaces in tool names
  if (!(ncol(pp) == 2)) {stop('Please ensure that there are no spaces in the tool names.')}
  
  names(pp) <- c("tool", "path")
  
  # check that all supplied paths exist
  pp$status <- file.exists(pp$path)
  
  if (all(pp$status)) { 
    message('All paths are valid')
    # convert all the tool names to lower case to avoid confusion
    pp <- plyr::mutate(pp, tool = tolower(tool))
    } else {
    message('Supplied file path does not exist.')
    message(sapply(pp %>% dplyr::filter(status == FALSE) %>% dplyr::select(path), paste, '\n'))
    message('Please check that supplied paths are correct and try again.')
    }
  
  # check that all the dependencies are executable in principle, i.e we are able to run -h or
  # process a small sample fasta file
  
  my_tools <- c('signalp2', 'signalp3', 'signalp4', 'targetp', 'tmhmm', 'wolfpsort')
  get_paths <- function(tool_name) {pp %>% filter(tool == tool_name) %>% select(path)}
  
  # siganlp2:
  
  if (suppressWarnings(system(paste(get_paths('signalp2'), '-h'), intern = TRUE)[2]) ==
      'Usage: signalp -t euk|gram+|gram- [options] seqfile') 
      { message('signalp2 test run completed.')
  } else {
        message('signalp2 test run failed; check if it is installed correctly.')
  }
  
  # signalp3:
  
  if (suppressWarnings(system(paste(get_paths('signalp3'), '-h'), intern = TRUE)[2]) ==
      'Usage: signalp -t euk|gram+|gram- [options] seqfile') 
      { message('signalp3 test run completed.')
  } else {
        message('signalp3 test run failed; check if it is installed correctly.')
  }
  
  # signalp4
  
  if (suppressWarnings(system(paste(get_paths('signalp4'), '-h'), intern = TRUE)[2]) ==
      '  Description: Predict signal peptide and cleavage site.') 
     { message('signalp4 test run completed.')
  } else {
       message('signalp4 test run failed; check if it is installed correctly.')
  }
  
  
 
  return(pp)
  
}



pa <- manage_paths(system.file("extdata", "sample_paths", package="SecretSanta"))
expect_match((system(paste((pa$path[1]), '-h'))), 'Description: Predict signal peptide and cleavage site.')
t1 <- system(paste((pa$path[1]), '-h'), intern = TRUE)[2]

#signalp4:

expect_match(t1[2], 'Description: Predict signal peptide and cleavage site.')

#signalp2:
signalp2

my_tools <- c('signalp2', 'signalp3', 'signalp4', 'targetp', 'tmhmm', 'wolfpsort')
get_paths <- function(tool_name) {pa %>% filter(tool == tool_name) %>% select(path)}

sapply(my_tools, get_paths)
  
pa %>% filter(tool == 'signalp2') %>% select(path)

system(paste(get_paths('signalp2'), '-h'))


