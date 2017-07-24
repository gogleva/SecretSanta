#' tmhmm function
#'
#' This function calls local TMHMM
#' Output:
#' "len=": the length of the protein sequence. \cr
#' "ExpAA=": The expected number of amino acids intransmembrane helices (see above). \cr
#' "First60=": The expected number of amino acids in transmembrane helices in the first 60 amino acids of the protein (see above). \cr
#' "PredHel=": The number of predicted transmembrane helices by N-best. \cr
#' "Topology=": The topology predicted by N-best. \cr
#' @param proteins input file with proteins
#' @param output_type what we want from the output
#' @export
#' @examples 
#' r1 <- tmhmm("SecretSanta/inst/extdata/sample_prot.fasta", 'test')

tmhmm <- function(proteins, output_type) {
  message("running TMHMM locally...")
  full_pa <- as.character(secret_paths %>% filter(tool == 'tmhmm') %>% select(path))
  result <- tibble::as.tibble(read.table(text = (system(paste(full_pa, proteins, '--short'), intern = TRUE))))
  names(result) <- c("gene_id", "length", "ExpAA",
                     "First60", "PredHel", "Topology")
  return(result)
}


