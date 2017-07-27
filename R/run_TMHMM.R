#' tmhmm function
#'
#' This function calls local TMHMM
#' Output: \cr
#' "len=": the length of the protein sequence. \cr
#' "ExpAA=": The expected number of amino acids intransmembrane helices (see above). \cr
#' "First60=": The expected number of amino acids in transmembrane helices in the first 60 amino acids of the protein (see above). \cr
#' "PredHel=": The number of predicted transmembrane helices by N-best. \cr
#' "Topology=": The topology predicted by N-best. \cr
#' @param input_obj input object, an instance of CBSResult class, prefered input should contain mature_fasta
#' @export
#' @examples 
#' r1 <- tmhmm("SecretSanta/inst/extdata/sample_prot.fasta", 'test')

# to do: make this function to produce an object of CBSResult class

tmhmm <- function(input_obj) {
  message("running TMHMM locally...")
  
  fasta <- getInfasta(input_obj) #placeholder, should be suitable for objects with mature_fasta slots only
  out_tmp <- tempfile()
  Biostrings::writeXStringSet(fasta, out_tmp)
  
  full_pa <- as.character(secret_paths %>% filter(tool == 'tmhmm') %>% select(path))
  result <- tibble::as.tibble(read.table(text = (system(paste(full_pa, out_tmp, '--short'), intern = TRUE))))
  names(result) <- c("gene_id", "length", "ExpAA",
                     "First60", "PredHel", "Topology")
  result <- (result %>% filter(PredHel == 'PredHel=0'))
  
  out_obj <- TMhmmResult(in_fasta = fasta,
                         out_fasta = fasta, #placeholder
                         tm_tibble = result)
  if (validObject(out_obj)) {return(out_obj)}
  }

### tests

#t1 <- tmhmm("SecretSanta/inst/extdata/sample_prot.fasta", 'test') # old version

t2 <- tmhmm(step1_sp2) # obj of SignalpResult class

