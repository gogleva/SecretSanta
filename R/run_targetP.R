#' tragetp function
#'
#' This function calls local targetp
#' Please check targetp README before running this fucntion:
#' TARGETP         'targetp-1.1' directory
#' TMP             location for temporary data
#' @param proteins input file with proteins
#' @param output_type what we want from the output
#' @export
#' @examples 

targetp <- function(proteins, output_type) {
  message("running targetp locally...")
  full_pa <- as.character(secret_paths %>% filter(tool == 'targetp') %>% select(path))
  result <- tibble::as.tibble(read.table(text = (system(paste(full_pa, "-N", proteins), intern = TRUE))))
}

#tests:

#tp <- targetp(proteins = "/home/anna/anna/Labjournal/SecretSanta/inst/extdata/sample_prot.fasta", output_type = 'something')
#system(paste(full_pa, "-N", "/home/anna/anna/Labjournal/SecretSanta/inst/extdata/sample_prot.fasta"), intern = TRUE)
#somethig wrong with targetp configuration, keeps returnig empty outputs





