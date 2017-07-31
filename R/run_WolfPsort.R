#' run_WolfPsort function
#'
#' This function runs WoLF PSORT to predict protein cellular sub-localisation
#' Recommended to run on the late stages of secretome prediction
#' @param input_obj Object of CSBResult class
#' @export
#' @examples

wolfpsort <- function(input_obj){
  if (is(input_obj, "CBSResult")) {} else {stop('input_object does not belong to CBSResult superclass')}
  if (length(getOutfasta(input_obj)) == 0) {stop('the input object contains empty out_fasta slot')}
  
  message("running WoLF PSORT locally...")
  
  fasta <- getOutfasta(input_obj)
  out_tmp <- tempfile()
  Biostrings::writeXStringSet(fasta, out_tmp)

  full_pa <- as.character(secret_paths %>% filter(tool == 'tmhmm') %>% select(path))
  
   # we assume/require that the object has filled outFasta slot already
  }
