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

m_slicer <- function(input_object, length_threshold) {
                    input_object <- aa
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
                      
                    # one AAStringSet object:
                    
                    many_slices <- function(n) {
                      sl <- unlist(sapply(X = unlist(smi[n]), FUN = slice, seq = input_object[n]))
                                               }
                    
                    smt <- sapply(1:length(smi), many_slices)
                    return(do.call(c, unlist(smt)))
                    
        }


l1 <- list(a = c(1:10), b = c(11:20))
l2 <- list(c = c(21:30), d = c(31:40))
# sum the corresponding elements of l1 and l2
mapply(sum, l1$a, l1$b, l2$c, l2$d)

listA <- input_object
listB <- smi

mapply(slice, listA[1], unlist(listB[1])
    
test_function <- function(n) {return(n + 10)}       
       
       
sapply(X = c(1:10), test_function)

       
