#' signalp function
#'
#' This function calls local SignalP. 
#' Please ensure that respective version of SignalP is downloaded and installed + added to $PATH variable
#' @param proteins input file/data frame with proteins
#' @param version  specify SignalP version to run; versions available: \cr
#'                 1.1 - http://www.cbs.dtu.dk/services/SignalP-1.1/ \cr
#'                 2.0 - http://www.cbs.dtu.dk/services/SignalP-2.0/ \cr
#'                 3.0 - http://www.cbs.dtu.dk/services/SignalP-3.0/ \cr
#'                 4.0 - http://www.cbs.dtu.dk/services/SignalP-4.0/ \cr
#'                 4.1 - http://www.cbs.dtu.dk/services/SignalP-4.1/ \cr
#' @export
#' @examples
#' signalp(proteins = "/SecretSanta/inst/extdata/sample_prot.fasta", version = 2)

signalp <- function(proteins, version) {
  message("running signalP locally...")
  result <- system(paste("signalp -t euk", proteins))
  
}


#' Parse_signalP_web function
#'
#' This function scraps web output of SignalP tool and returns output in short format organised in a dataframe
#' @param url url used to produce the output
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




