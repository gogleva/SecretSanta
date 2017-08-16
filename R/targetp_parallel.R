#combine_TaretpResult function
#'
#' This helper function combines multiple instances of TargetpResult class, typically generated with parLapply
#' @param arguments - a list of TargetpResult objects to be combined in one
#' @export
#' @examples 
#' inp2 <- CBSResult(in_fasta = readAAStringSet(system.file("extdata", "tail_prot.fasta", package = "SecretSanta")))
#' inp3 <- CBSResult(in_fasta = readAAStringSet(system.file("extdata", "tail2_prot.fasta", package = "SecretSanta")))
#' inp4 <- CBSResult(in_fasta = readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta")))
#' 
#' tp1 <- targetp(input_obj = inp2, network_type = 'N', run_mode = 'starter', paths = my_pa)
#' tp2 <- targetp(input_obj = inp3, network_type = 'N', run_mode = 'starter', paths = my_pa)
#' tp3 <- targetp(input_obj = inp4, network_type = 'N', run_mode = 'starter', paths = my_pa)
#' 
#' obj <- list(tp1, tp2, tp3)
#  combined_tp <- combine_TargetpResult(obj)


combine_TargetpResult <- function(arguments) {
  if (all(sapply(arguments, is, 'TargetpResult'))) {
  } else {                               
    stop('Some objects from arguments list do not belong to TargetpResult class.')
  }
  
  c_in_fasta <- do.call(c, (lapply(arguments, getInfasta)))
  c_out_fasta <- do.call(c, (lapply(arguments, getOutfasta)))
  c_tp_tibble <- do.call(rbind, (lapply(arguments, getTPtibble)))
  
  c_obj <- TargetpResult(in_fasta = c_in_fasta,
                         out_fasta = c_out_fasta,
                         tp_tibble = c_tp_tibble)
                         }

#' targetp function
#'
#' This function calls local targetp to predict subcellular localisation of a protein.
#' Parallelised
#' @param input_object    an instance of CBSResult class containing protein sequences as on of the attributes
#' @param network_type possible values: 
#' \itemize{
#' \item P - for plants;
#' \item N - for non-plants;
#' }
#' @param run_mode
#' \itemize{
#' \item starter - if it is the first step in pipeline;
#' \item piper - if you run this function on the output of other CBS tools;
#' }
#' @param paths   tibble with paths to external dependencies, generated with \code{\link{manage_paths}} function
#' @return an object of TargetpResult class
#' @export
#' my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), use.names = TRUE)
# 
#' # assign this object to the input_fasta slot of SignalpResult object
#' inp <- CBSResult(in_fasta = aa)
#' 
#' tp_result <- targetp(input_object = inp, network_type = 'N', run_mode = 'starter', paths = my_pa)


targetp <- function(input_object, network_type, run_mode, paths) {
  
  # helper function: crop long names for AAStringSet object, return character vector
  crop_names <- function(x){unlist(stringr::str_split(x, " "))[1]}
  
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
  
  # checked that specified networks are valid  
  allowed_networks = c('P', 'N')
  
  if (network_type %in% allowed_networks) {
    message("running targetp locally...")
  } else {
    stop('Specified network_type is invalid.')  
  }
  
# simple function to run targetp on relatively small input files ~1K proteins
  
  simple_targetp <- function(aaSet){
      #----- Run targetp prediction:
      message(paste('Number of submitted sequences...', length(aaSet)))
  
      # convert fasta to a temporary file:
      out_tmp <- tempfile() #create a temporary file for fasta
      Biostrings::writeXStringSet(aaSet, out_tmp) #write tmp fasta file
  
      # get path to targetp executable
      full_pa <- as.character(paths %>% dplyr::filter(tool == 'targetp') %>% dplyr::select(path))
  
      # prep fasta:
      # generate cropped names for input fasta
      cropped_names <- unname(sapply(names(aaSet), crop_names))
      # replace long names with cropped names
      names(aaSet) <- cropped_names
  
      # prep networks argument:
      NN <- paste('-', network_type, sep = '')
  
      #run targetp:
      tp <- tibble::as.tibble(read.table(text = (system(paste(full_pa, NN, out_tmp), intern = TRUE)[1: length(aaSet) + 8])))
  
      if (network_type == 'N') {
        names(tp) <- c('gene_id', 'length', 'mTP', 'sp', 'other', 'TP_localization', 'RC')}
      else if (network_type == 'P') {
        names(tp) <- c('gene_id', 'length', 'cTP', 'mTP', 'sp', 'other', 'TP_localization', 'RC')  
      }
  
      tp <- tp %>% dplyr::filter(TP_localization == 'S')
      message(paste('Number of candidate secreted sequences', nrow(tp)))         
  
      candidate_ids <- tp %>% dplyr::select(gene_id) %>% unlist(use.names = FALSE)
      out_fasta_tp <- aaSet[candidate_ids]
  
      # generate output object:
      out_obj <- TargetpResult(in_fasta = aaSet,
                               out_fasta = out_fasta_tp,
                               tp_tibble = tp)
  
      if (validObject(out_obj)) {return(out_obj)}
  }
  

   # # Check input file size and decide how to run targetp: in parallel mode or not:
   if (length(fasta) <= 1000) {
     message('Ok for single processing')
     return(simple_targetp(fasta))
   } else {
     message('Input fasta contains >1000 sequences, entering batch mode...')
     message(paste('Number of submitted sequences...', length(fasta)))
     # split fasta:
     split_fasta <- split_XStringSet(fasta, 1000)
     
     # initiate cluster on all the cores available and supply required environment:
     no_cores <- detectCores()
     cl <- makeCluster(no_cores)
     clusterEvalQ(cl, library("SecretSanta"))
     
     clusterExport(cl=cl, varlist=c("my_pa"), envir=environment()) #  works oly with my_pa, not paths!
     
     # produces error: Error in get(name, envir = envir) : object 'signalp4' not found
     
     # run parallel targetp:
     result <- parLapply(cl, split_fasta, simple_targetp)
     stopCluster(cl)
     
     res_comb <- do.call(c,result)
     combined_TargetpResult <- combine_TargetpResult(unname(res_comb))
     
     # may be add some wraings/messages here:
     
     tp_count <- nrow(getTPtibble(combined_TargetpResult))    
     message(paste('Number of candidate secerted sequences...', tp_count))
     if (tp_count == 0) {warning('Targetp prediction yeilded 0 extracellular candidates')}
     
     closeAllConnections()
     return(combined_TargetpResult)
     
   } 
     
}


