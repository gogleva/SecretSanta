#' predict transmembrane domains with TOPCONS
#'
#' This function calls local TOPCONS to predict transmembarne domains in protein
#' sequence. \cr
#' \cr
#' @param input_obj an instance of SignalpResult class, \cr
#'                  input should contain mature_fasta slot;
#' @param TM  allowed number of TM domains in mature peptides,
#' recommended value <= 1; use TM = 0 for strict filtering                   
#' @param paths if topcons is not acessible globally, a file
#' conatining a full path to it's executable should be provided; for details
#' please check SecretSanta vignette.
#' @export
#' @return TopconsResult object
#' @examples 

topcons <- function(){}

