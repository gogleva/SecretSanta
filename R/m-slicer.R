#' m_slicer function
#'
#' Experimental option.
#' This function generates all possible subsequences starting with M.
#' Assumption: translation start sites might be mis-predicted in the original set of proteins.
#' Output of this step can be used as an input for secretome prediction pipeline to rescue secreted proteins
#' with mis-predicted start sites. 
#' 
#' @param input_object    an instance of CBSResult class or AAStringSet class containing protein sequences as on of the     attributes
#' @param len_threshold   sliced sequences below this threshold will be discarded
#' @param run_mode
#' \itemize {
#'   \item  \strong{slice} - to just slice input fasta, regardless of it's origin;
#'   \item  \strong{rescue} - to get proteins not predicted to be secreted on the initial run, generate slices; 
#'   }
#' @export
#' @examples 
#' # Example 1: generate proteins with alterative translation start site for AAStringSet object
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), use.names = TRUE)
#' m_slicer(aa, 100, run_mode = 'slice')
#' 
#' # Example 2: generate proteins with alterative translation start site for CBSResult object
#' my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
#' inp <- SignalpResult()
#' inp <- setInfasta(inp, aa)
#' s1_sp2 <- signalp(inp, version = 2, 'euk', run_mode = "starter", paths = my_pa)
#' slices <- m_slicer(s1_sp2, length_threshold = 100, run_mode = 'rescue')
#' inp_slices <- CBSResult(in_fasta = slices)
#' s2_sp2_rescue <- signalp(inp_slices, version = 2, 'euk', run_mode = 'starter', paths = my_pa)# => throws an error! too many sequences


m_slicer <- function(input_object, length_threshold, run_mode) {
  
                    # helper functions:
                    '%!in%' <- function(x,y)!('%in%'(x,y))
                    
                    # check that inputs are valid
                    
                    if (is(input_object, 'AAStringSet')) {
                      input_object <- input_object
                      if (run_mode != 'slice') {
                        stop("Please use run_mode 'slice' for an input object of AAStringSet class")
                        } 
                    } else if (is(input_object, 'CBSResult')) {
                      if (run_mode != 'rescue') {
                        stop("Please use run_mode 'rescue' for an input object of CBSResult class")
                      }
                      infa <- getInfasta(input_object)
                      outfa <- getOutfasta(input_object)
                      input_object <- infa[infa %!in% outfa]
                    }
  
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

'%!in%' <- function(x,y)!('%in%'(x,y))
