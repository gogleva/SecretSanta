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
  # ----- Check that inputs are valid
  
  # check that input object belong to CBSResult class
  
  if (is(input_obj, "CBSResult")) {} else {stop('input_object does not belong to CBSResult superclass')}
  
  # check that supplied runnig mode is valid
  
  if (run_mode %in% c('piper', 'starter')) {} else {stop("Run mode is invalid. Please use 'starter' to initiate prediction pipelie or 'piper' to continue")}
  
  # check that input_object contains non-empty in/out_fasta for starter/piper
  
  if (run_mode == 'starter') {
    if (length(getInfasta(input_obj)) != 0) {
      fasta <- getInfasta(input_obj)
    } else {stop('in_fasta attribute is empty')}
  } else if (run_mode == 'piper') {
    if (length(getOutfasta(input_obj)) != 0) {
      fasta <- getOutfasta(input_obj)
    } else {stop('out_fasta attribute is empty')}
  }
  
  allowed_networks = c('P', 'N')
  
  if (network_type %in% allowed_networks) {
    message("running targetp locally...")
  } else {
    stop('Specified network_type is invalid.')  
  }
    
  #----- Run targetp prediction:
  
  full_pa <- as.character(paths %>% filter(tool == 'targetp') %>% select(path))
  result <- tibble::as.tibble(read.table(text = (system(paste(full_pa, "-N", proteins), intern = TRUE))))
}

#tests:

#tp <- targetp(proteins = "/home/anna/anna/Labjournal/SecretSanta/inst/extdata/sample_prot.fasta", output_type = 'something')
#system(paste(full_pa, "-N", "/home/anna/anna/Labjournal/SecretSanta/inst/extdata/sample_prot.fasta"), intern = TRUE)
#somethig wrong with targetp configuration, keeps returnig empty outputs





