#' manage_paths function
#'
#' This function calls manages pathways for signalP versions
#' @param path_file  file with listed paths for different signalp versions
#' @export
#' @examples
#' 

manage_paths <- function(path_file) {
  # read path file in a tibble
  pp <- readr::read_delim(path_file, delim = ' ', col_names = FALSE)
  names(pp) <- c("tool", "path")
  # now check that all supplied paths exist
  pp$status <- file.exists(pp$path)
  if (any(pp$status == FALSE)) { 
    message('supplied file path does not exist')} else {
      return(pp)
    }
}

testmm


