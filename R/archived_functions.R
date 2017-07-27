#' Parse_signalP_web function
#'
#' This function scraps web output of SignalP tool and returns output in short format organised in a dataframe
#' @param url url used to produce the output
#' @export
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

#' pipe_signalp function
#'
#' This function pipes outputs from different version of signalp
#' @param piping_sequence  specify piping sequence: vector of signalp versions \cr
#'                         accepted formats: c(1,2,3) or c("signalp2", "signalp3", "signalp4")
#' @param input_fasta fasta file with proteins to initiate the pipe
#' @examples
#' piping_seq <- c("signalp2", "signalp3", "signalp4")
#' piping_seq_num <- c(2,3,3)
#' pipe_signalp(piping_sequence = piping_seq, 'some_fasta')
#' pipe_signalp(piping_sequence = piping_seq_num, 'some_fasta')
#'

pipe_signalp <- function(piping_sequence, input_fasta) {
  # check that pipong sequence does not contain duplicates
  if (length(unique(piping_sequence)) == length(piping_sequence)) {
    # generate starting message
    message("Hi, let's pipe!")
    if (all(stringr::str_detect(piping_sequence, "signalp"))) {
      piping_sequence
    } else if(is.numeric(piping_sequence)){
      piping_sequence <- paste('signalp', piping_seq_num, sep = '')
    }
    message(paste(piping_sequence[-length(piping_sequence)], '--> '), piping_sequence[length(piping_sequence)])
    message('pipeline includes ', length(piping_seq), ' steps')
  }else{
    message('Aborting..')
    stop('please make sure that there are no duplicte tools in the specified piping sequence')
  }  
  
  # piper:
  
  # run signalp#1(input_fasta) -> outputs tibble -> generate (fasta)'
  # run signalp#2(fasta') -> output tibble' -> generate (fasta)''
  # run signalp#3(fasta'') -> output tibble'' == result
  
  
}

#helper recursive piper function: 

#piper <- function(some_sequence) {
#   if  
#}


#to do: teach run_signalp to output fasta file and/or tibble with results

#' Convert \code{data.frame} to \code{list}.
#' 
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @param x A \code{data.frame} object.
#' @examples
#' my_result <- foo(iris)
#'
foo <- function(x) {
  x %>%
    as.list()
}


