#' check_khdel function
#'
#' This function checks presence of terminal KDEL/HDEL sequences in the candidate secreted proteins
#' @param input_obj input object of CBSResult superclass
#' @param run_mode 'starter' or 'piper' 
#' @export
#' @examples 
#' check_khdel(step1_sp2, run_mode = 'starter')

check_khdel <- function(input_obj, run_mode) {
  # check the input
  if (is(input_obj, "CBSResult")) {} else stop('input_object does not belong to CBSResult superclass')
  
  # starting message:
  message("checking for terminal ER retention signals...")
  
  # determine which fasta to take
  if (run_mode == 'piper') fasta <- getOutfasta(input_obj) else fasta <- getInfasta(input_obj)
  
  # find and remove termial ER retention motifs
  ER1 <- AAString('KDEL')
  ER2 <- AAString('HDEL')
  tails <- subseq(fasta, -4, -1) # last 4 aminno acids in each protein
  un <- !(as.logical(vcountPattern(ER1, tails)) | as.logical(vcountPattern(ER2, tails)))
  
  # fasta without terminal KDELs/HDEls
  non_retained <- fasta[un] 
  
  out_obj <- ErResult(in_fasta = fasta,
                      out_fasta = non_retained,
                      retained = fasta[!un])
  
  if (validObject(out_obj)) {return(out_obj)}
}



