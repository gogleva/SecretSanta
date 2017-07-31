#' run_WolfPsort function
#'
#' This function runs WoLF PSORT to predict protein cellular sub-localisation
#' And returns the most probbale
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
  full_pa <- as.character(secret_paths %>% filter(tool == 'wolfpsort') %>% select(path))
  wolf <- system(paste(full_pa, 'fungi <', out_tmp), intern = TRUE)
  
  #parse wolf output
  clean_strings <- function(x, field){unlist((stringr::str_split(x, " ")))[c(field)]}
  
  gene_id <- sapply(X = wolf, field = 1, clean_strings, USE.NAMES = FALSE)
  localization <- sapply(X = wolf, field = 2, clean_strings, USE.NAMES = FALSE)
  
  #assemble result tibble with gene id and most probable sibsellular localisation
  wolf_tbl <- as_tibble(data.frame(gene_id, localization)) %>% filter(gene_id != '#') %>% filter(localization == 'extr')
  
  #assemble wolf result object:
  out_obj <- WolfResult(in_fasta = fasta,
                        out_fasta = fasta[wolf_tbl$gene_id], #place holder
                        wolf_tibble = wolf_tbl)
  
  if (validObject(out_obj)) {return(out_obj)}
}


#tests:

fafa <- step1_sp2
w <- wolfpsort(fafa)


