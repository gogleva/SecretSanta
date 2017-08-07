#' signalp function
#'
#' This function calls local targetp to predict subcellular localisation of a protein.
#' @param input_object    an instance of CBSResult class containing protein sequences as on of the attributes
#' @param network_type possible values: 
#' \itemize{
#' \item P - for plants;
#' \item N - for non-plants;
#' @param run_mode
#' \itemize{
#' \item starter - if it is the first step in pipeline;
#' \item piper - if you run this function on the output of other CBS tools;
#' }
#' @param paths   tibble with paths to external dependencies, generated with \code{\link{manage_paths}} function
#' @return an object of SignalpResult class
#' @export

targetp <- function(input_object, network_type, run_mode, paths) {
  message("running targetp locally...")
  full_pa <- as.character(paths %>% filter(tool == 'targetp') %>% select(path))
  result <- tibble::as.tibble(read.table(text = (system(paste(full_pa, "-N", proteins), intern = TRUE))))
}

#tests:

#tp <- targetp(proteins = "/home/anna/anna/Labjournal/SecretSanta/inst/extdata/sample_prot.fasta", output_type = 'something')
#system(paste(full_pa, "-N", "/home/anna/anna/Labjournal/SecretSanta/inst/extdata/sample_prot.fasta"), intern = TRUE)
#somethig wrong with targetp configuration, keeps returnig empty outputs





