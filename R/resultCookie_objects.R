#' An S4 class to represent intermediate and final outputs of the
#' secretome prediction pipeline
#' 
#' @slot fasta initial fasta file
#' @slot mature_fasta fasta with mature sequences
#' 
#' 
#' 
resultCookie <- setClass("resultCookie",
                    slots = list(fasta = "character",
                                 mature_fasta = "character"))
