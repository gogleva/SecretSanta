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
  
  return(pp)
}


#alias signalp4="/home/anna/anna/Labjournal/SecretSanta_external/signalp-4.1/signalp"
#alias signalp3="/home/anna/anna/Labjournal/SecretSanta_external/signalp-3.0/signalp"
#alias signalp2="/home/anna/anna/Labjournal/SecretSanta_external/signalp-2.0/signalp"

testmm

my_check <- function(x) file.exists(x) 


if file.exists(p) {}
