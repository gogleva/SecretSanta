#' check_khdel function
#'
#' This function checks for terminal KDEL/HDEL sequences in the candidate secreted proteins
#' @param proteins input file with proteins
#' @param input_type path, object or something else?
#' @param output_type what we want from the output: fasta, tibble? list of ids?
#' @export
#' @examples 
#' check_khdel("/home/anna/anna/Labjournal/SecretSanta/inst/extdata/sample_prot.fasta", 'some')

check_khdel <- function(proteins, output_type) {
  message("checking for terminal ER retention signals...")
  ff <- read.fasta(proteins)
  # helper function to gen n tail bases
  tail_bases <- function(x){
    L <- length(x)
    return(paste(x[(L - 3):L], collapse = ''))
  }

  if (any(which(mapply(tail_bases, ff) == 'hdel')) | any(which(mapply(tail_bases, ff) == 'kdel'))) {
    message('oh no, kdel')
    res <- attr(c(which(mapply(tail_bases, ff) == 'hdel'), which(mapply(tail_bases, ff) == 'kdel')), "name")
    print(res)
    }else{
      message('all good!')
    }
}


