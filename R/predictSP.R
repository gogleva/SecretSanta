#' PredictSP function
#'
#' This function call SignalP for prediction of signal peptides for amino acid sequences
#' @param proteins input file/data frame with proteins
#' @param version  specify SignalP version to run
#' @param mode  SignalP running mode, web or standalone (require local installation)
#' @keywords signal peptide, proteins, signal peptide prediction, HMM
#' @export
#' @examples
#' predictSP()

predictSP <- function(proteins, version, mode) {
  print(paste("hi", proteins, "let's get started!", sep = " "))
}