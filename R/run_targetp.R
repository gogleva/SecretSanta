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
  
  if (is(input_object, "CBSResult")) {} else {stop('input_object does not belong to CBSResult superclass')}
  
  # check that supplied runnig mode is valid
  
  if (run_mode %in% c('piper', 'starter')) {} else {stop("Run mode is invalid. Please use 'starter' to initiate prediction pipelie or 'piper' to continue")}
  
  # check that input_object contains non-empty in/out_fasta for starter/piper
  
  if (run_mode == 'starter') {
    if (length(getInfasta(input_object)) != 0) {
      fasta <- getInfasta(input_object)
    } else {stop('in_fasta attribute is empty')}
  } else if (run_mode == 'piper') {
    if (length(getOutfasta(input_object)) != 0) {
      fasta <- getOutfasta(input_object)
    } else {stop('out_fasta attribute is empty')}
  }
  
  allowed_networks = c('P', 'N')
  
  if (network_type %in% allowed_networks) {
    message("running targetp locally...")
  } else {
    stop('Specified network_type is invalid.')  
  }
    
  #----- Run targetp prediction:

  # convert fasta to a temporary file:
  out_tmp <- tempfile() #create a temporary file for fasta
  Biostrings::writeXStringSet(fasta, out_tmp) #write tmp fasta file
  
  # make a system call of targetp based on the tmp file
  
  full_pa <- as.character(paths %>% dplyr::filter(tool == 'targetp') %>% dplyr::select(path))
  
  # helper function: crop long names for AAStringSet object, return character vector
  crop_names <- function(x){unlist(stringr::str_split(x, " "))[1]}
  message(paste('Number of submitted sequences...', length(fasta)))
  
  # generate cropped names for input fasta
  cropped_names <- unname(sapply(names(fasta), crop_names))
  # replace long names with cropped names
  names(fasta) <- cropped_names
  
  #run targetp:
  
  NN <- paste('-', network_type, sep = '')
  tp <- tibble::as.tibble(read.table(text = (system(paste(full_pa, NN, out_tmp), intern = TRUE)[1: length(fasta) + 8])))
  names(tp) <- c('gene_id', 'length', 'mTP', 'sp', 'other', 'TP_localization', 'RC')
  tp <- tp %>% dplyr::filter(TP_localization == 'S')
  message(paste('Number of candidate secreted sequences', nrow(tp)))         
  
  # generate output object:
  
#  out_obj <- TargetpResult()

  return(tp)
}


# tests:
# 
# my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
# 
# # initialise SignalpResult object
# inp <- SignalpResult()
# 
# # read fasta file in AAStringSet object
# aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), use.names = TRUE)
# 
# # assign this object to the input_fasta slot of SignalpResult object
# inp <- setInfasta(inp, aa)
# 
# 
# test <- targetp(input_object = inp, network_type = 'N', run_mode = 'starter', paths = my_pa)








