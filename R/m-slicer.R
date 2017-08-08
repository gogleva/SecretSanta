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
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), use.names = TRUE)
#' 

m_slicer <- function(input_object, length_threshold, n) {
                    mi <- vmatchPattern('M', input_object)
                    smi <- startIndex(mi) #all M-positions
                    fifi <- function(x) { unlist(x)[unlist(x) >1]}
                    
                    # remove 1's from smi:
                    smi <- lapply(smi, fifi)
                    wi <- which(lengths(smi) > 0)
                    input_object <- input_object[wi]
                    smi <- smi[lengths(smi) > 0]
  
                    # slice one AAString
                    slice <- function(x, seq) {
                      st <- subseq(seq, start = x, end = -1)
                      names(st) <- paste(unlist(strsplit(names(st), ' '))[1], '_slice_M', x, sep = '')
                      if (width(st) >= length_threshold) {return (st)}
                                }
                                              }
  
                    # one AAStringSet object:
                    #slices <- do.call(c, unlist(sapply(X = unlist(smi[n]), FUN = slice, seq = input_object[n])))
                    #return(slices)

                    listA <- input_object
                    listB <- smi
                    
                    mapply(function(X,Y) {
                                          sapply(1:length(X), function(x) )
                    })
                    
                    many_slices <- function(x,y) {
                      result <- unlist(sapply(X = unlist(x), FUN = slice, seq = y))
                      return(result)
                    }
                    
                    
                    listA <- list(matrix(rnorm(2000), nrow=10),
                                  matrix(rnorm(2000), nrow=10))
                    listB <- list(matrix(rnorm(2000), nrow=10),
                                  matrix(rnorm(2000), nrow=10))
                    mapply(function(X,Y) {
                      sapply(1:10, function(row) cor(X[row,], Y[row,]))
                    }, X=listA, Y=listB)

  
        }
