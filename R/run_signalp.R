#' signalp function
#'
#' This function calls local SignalP
#' @param proteins input file/data frame with proteins
#' @param version  specify SignalP version to run
#' @keywords signalp
#' @export
#' @examples
#' signalp(proteins = "/SecretSanta/inst/extdata/sample_prot.fasta", mode = 'local', version = 2)

signalp <- function(proteins, version) {
  print("running signalP locally")  
  print(paste("/home/anna/anna/soft/signalp_2.0/signalp-2.0/signalp -t euk", proteins))
}


#' Parse_signalP_web function
#'
#' This function scraps web output of SignalP tool and returns output in short format organised in a dataframe
#' @param url url used to produce the output
#' @keywords signal peptide, proteins, signal peptide prediction, HMM, web scrapping
#' @export
#' @examples
#' Parse_signalP_web()

Parse_signalP_web <- function(url){
  webpage <- xml2::read_html(url)
  raw_data_html <- rvest::html_nodes(webpage, 'pre')
  data <- rvest::html_text(raw_data_html)
  raw <- unlist(strsplit(data, '\n'))
  le <- length(raw)
  d <- raw[9:le - 4] #draft parameters
  d2 <- as.data.frame(do.call("rbind", d2)) 
  d2
}




