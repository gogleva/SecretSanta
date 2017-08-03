#' check_khdel function
#'
#' This function checks presence of terminal KDEL/HDEL sequences in the candidate secreted proteins.
#' @param input_obj input object of CBSResult class
#' @param run_mode 
#' \itemize{
#' \item starter - if it is the first step in pipeline;
#' \item piper - if you run this function on the output of other CBS tools;
#' }
#' @export
#' @examples 
#' # check ER retention signals in CBSResult object before running signalp or any other predictions
#' inp <- SignalpResult()
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), use.names = TRUE)
#' inp <- setInfasta(inp, aa)
#' et_s <- check_khdel(inp, run_mode = 'starter')
#' 
#' # check ER retention signal in the signalp output, 'starter' mode
#' et_sp <- check_khdel(step1_sp2, run_mode = 'starter')
#' 
#' # check ER retention signal in the signalp output, 'piper' mode
#' et_piper <- check_khdel(step1_sp2, run_mode = 'piper')

check_khdel <- function(input_obj, run_mode) {
  # check the input
  if (is(input_obj, "CBSResult")) {} else stop('input_object does not belong to CBSResult superclass')
  
  # starting message:
  message("checking for terminal ER retention signals...")
  
  # determine which fasta to take
  if (run_mode == 'piper') fasta <- getOutfasta(input_obj) else fasta <- getInfasta(input_obj)
  
  message(paste('Number of submitted sequences...', length(fasta)))
  
  if (length(fasta) == 0) {stop('query fasta is empty, please ensure you are using correct run_mode')}
  
  # find and remove termial ER retention motifs
  ER1 <- AAString('KDEL')
  ER2 <- AAString('HDEL')
  tails <- subseq(fasta, -4, -1) # last 4 aminno acids in each protein
  un <- !(as.logical(vcountPattern(ER1, tails)) | as.logical(vcountPattern(ER2, tails)))
  
  # fasta without terminal KDELs/HDELs
  non_retained <- fasta[un] 
  
  out_obj <- ErResult(in_fasta = fasta,
                      out_fasta = non_retained,
                      retained = fasta[!un])
  
  ret_count <- length(fasta[!un])
  message(paste('Number of sequences with terminal ER retention signals detected...', ret_count))
  
  
  if (validObject(out_obj)) {return(out_obj)}
}



