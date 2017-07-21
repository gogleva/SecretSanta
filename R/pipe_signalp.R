#' pipe_signalp function
#'
#' This function pipes outputs from different version of signalp
#' @param piping_sequence  specify piping sequence: vector of signalp versions \cr
#'                         accepted formats: c(1,2,3) or c("signalp2", "signalp3", "signalp4")
#' @param input_fasta fasta file with proteins to initiate the pipe
#' @export
#' @examples
#' piping_seq <- c("signalp2", "signalp3", "signalp4")
#' piping_seq_num <- c(2,3,3)
#' pipe_signalp(piping_sequence = piping_seq, 'some_fasta')
#' pipe_signalp(piping_sequence = piping_seq_num, 'some_fasta')
#'

pipe_signalp <- function(piping_sequence, input_fasta) {
  # generate starting message
  message("Hi, let's pipe!")
  if (all(stringr::str_detect(piping_sequence, "signalp"))) {
    piping_sequence
  } else if(is.numeric(piping_sequence)){
    piping_sequence <- paste('signalp', piping_seq_num, sep = '')
  }
  message(paste(piping_sequence[-length(piping_sequence)], '--> '), piping_sequence[length(piping_sequence)])
  # do the piping
}



