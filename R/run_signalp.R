#' decideSP function
#'
#' This function gathers output and decides witch signalP wrapper to run
#' @param version  specify SignalP version to run
#' @param mode  SignalP running mode, web or standalone (require local installation)
#' @keywords signal peptide, proteins, signal peptide prediction, HMM
#' @export
#' @examples
#' decideSP()


decideSP <- function(version, mode) {
  pass
}

#' PredictSP function
#'
#' This function calls local SignalP
#' @param proteins input file/data frame with proteins
#' @param version  specify SignalP version to run
#' @param mode  SignalP running mode, web or standalone (require local installation)
#' @keywords signal peptide, proteins, signal peptide prediction, HMM
#' @export
#' @examples
#' predictSP()


predictSP <- function(proteins, version, mode) {
  print(paste("hi", proteins, "let's get started!", sep = " "))
  system("/home/anna/anna/soft/signalp_2.0/signalp-2.0/signalp -h")
}

#' Submit_signalP_web function
#'
#' This function submits input to SignalP web tool and returns the actual url which could be used by Parse_signalP_web
#' @param fasta fasta file with protein sequences
#' @param mode  specify output format
#' @keywords signal peptide, proteins, signal peptide prediction, HMM, web scrapping
#' @export
#' @examples
#' Submit_signalP_web()
 
Submit_signalP_web <- function(fasta, mode) {
  pass
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




