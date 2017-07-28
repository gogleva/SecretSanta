#' check_khdel function
#'
#' This function checks presence of terminal KDEL/HDEL sequences in the candidate secreted proteins
#' @param input_obj input object of CBSResult superclass
#' @export run_mode 'starter' or 'piper' 
#' @examples 
#' check_khdel("/home/anna/anna/Labjournal/SecretSanta/inst/extdata/sample_prot.fasta", 'some')

check_khdel <- function(input_obj, run_mode) {
  
  # check the input
  if (is(input_obj, "CBSResult")) {} else stop('input_object does not belong to CBSResult superclass')
  
  # starting message:
  message("checking for terminal ER retention signals...")
  
  # determine which fasta to take
  if (run_mode == 'piper') fasta <- getOutfasta(input_obj) else fasta <- getInfasta(input_obj)
  
  ER1 <- AAString('KDEL')
  ER2 <- AAString('HDEL')
  
  tails <- subseq(fasta, -4, -1) # last 4 aminno acids in each protein
  
  un <- !(as.logical(vcountPattern(ER1, tails)) | as.logical(vcountPattern(ER2, tails)))
  #un <- !(qu1 | qu2)
  
  return(fasta[un])
}


## Tests:
check_khdel(step1_sp2, run_mode = 'starter')
