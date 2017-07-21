#' manage_paths function
#'
#' This function calls manages pathways for CBS tools: signalp, targetp ... etc
#' @param path_file  2-column space-separated file with listed paths for CBS tools
#' @export
#' @examples
#' secret_paths <- manage_paths("SecretSanta/inst/extdata/sample_paths")

manage_paths <- function(path_file) {
  # read path file in a tibble
  pp <- readr::read_delim(path_file, delim = ' ', col_names = FALSE)
  names(pp) <- c("tool", "path")
  # now check that all supplied paths exist
  pp$status <- file.exists(pp$path)
  if (all(pp$status)) { 
    return(pp)
    } else {
    message('Warning! supplied file path does not exist')
    message(sapply(pp %>% filter(status == FALSE) %>% select(path), paste, '\n'))
    message('Please check the paths and try again')
    }
}