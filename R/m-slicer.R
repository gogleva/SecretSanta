#' m_slicer function
#'
#' Experimental option
#' This function generates all possible subsequences starting with M.
#' Assumption: translation start sites might be mis-predicted in the original set of proteins.
#' Output of this step can be used as an iput for secretome prediction pipeline
#' 
#' @param input_object    an instance of CBSResult class or AAStringSet class containing protein sequences as on of the     attributes
#' @param len_threshold   sliced sequences below this threshold will be discarded
#' @export
#' @examples 


m_slicer <- function(input_object, length_threshold) {
                    x <- 1      
                          }

# sample AAstringSet:

aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), use.names = TRUE)

mi <- vmatchPattern('M', aa)
smi <- startIndex(mi) #all M-positions

slice <- function(x, seq) {subseq(seq, start = x, end = -1)}

sapply(X = unlist(smi[3]), FUN = slice, seq = aa[3])




