#' m_slicer function
#'
#' NB: this is an \strong{experimental option}.\cr
#' \cr
#' This function generates all possible subsequences (slices) starting with 
#' methionines (M). This might be usefull when translation start sites assumed
#' to be mis-predicted for some of the provided proteins. For example, in cases
#' when the set is obtained after \emph{de novo} genome or transcriptome
#' assembly.
#' \cr
#' \cr
#' Output of this step can be used as an input for secretome prediction
#' pipeline to rescue secreted proteins with potentially mis-predicted start
#' sites. Please proceed with caution.
#' 
#' @param input_obj        an instance of CBSResult class or AAStringSet
#'    class containing protein sequences as on of the attributes
#' @param min_len     sliced sequences below this threshold will be
#'    discarded
#' @param run_mode    \strong{slice} - to just slice input fasta, regardless
#' of its origin; \cr
#' \strong{rescue} - to get proteins not predicted to be secreted on the 
#' initial run, generate slices; 
#' @return a set of sliced sequences, AAStringSet object
#' @export     
#' @examples 
#' # Example 1: generate proteins with alterative translation start site for
#' # AAStringSet object
#' aa <- readAAStringSet(system.file("extdata","sample_prot_100.fasta",
#' package = "SecretSanta"))
#' m_slicer(aa[1:10], 100, run_mode = 'slice')
#' 
#' # Example 2: generate proteins with alterative translation start site for
#' # CBSResult object
#' inp <- CBSResult(in_fasta = aa[1:10])
#' s1_sp2 <- signalp(inp, version = 2, organism = 'euk',
#' run_mode = "starter")
#' slices <- m_slicer(s1_sp2, min_len = 100, run_mode = 'rescue')
#' inp_slices <- CBSResult(in_fasta = slices)
#' s2_sp2_rescue <- signalp(inp_slices, version = 2, organism = 'euk',
#' run_mode = 'starter')

m_slicer <- function(input_obj, min_len, run_mode = c('slice', 'rescue')) {
    
    # helper function:
    '%!in%' <- function(x, y) !('%in%'(x, y))
    
    # check that inputs are valid
    if (missing(run_mode)) {stop('missing argument: run_mode')}
    run_mode = match.arg(run_mode)
    
    
    if (is(input_obj, 'AAStringSet') &&
            (run_mode != 'slice')) {
        stop("Use run_mode 'slice' for an input object of AAStringSet class")
    }
    
    if (is(input_obj, 'CBSResult') &&
            (run_mode != 'rescue'))
    {
        stop("Use run_mode 'rescue' for an input object of CBSResult class")
    }
    
    # hadle inputs
    
    if (is(input_obj, 'AAStringSet')) {
        input_obj <- input_obj
    }
    
    if (is(input_obj, 'CBSResult')) {
        infa <- getInfasta(input_obj)
        outfa <- getOutfasta(input_obj)
        input_obj <-
            infa[names(infa) %!in% names(outfa)]
    }
    
    mi <- vmatchPattern('M', input_obj)
    smi <- startIndex(mi) #all M-positions
    fifi <- function(x) {unlist(x)[unlist(x) > 1]}
    
    # remove 1's from smi:
    smi <- lapply(smi, fifi)
    wi <- which(lengths(smi) > 0)
    input_obj <- input_obj[wi]
    smi <- smi[lengths(smi) > 0]
    
    # slice one AAString
    slice <- function(x, seq) {
        st <- subseq(seq, start = x, end = -1)
        names(st) <-
            paste(unlist(strsplit(names(st), ' '))[1],
                        '_slice_M',
                        x,
                        sep = '')
        if (width(st) >= min_len) {
            return (st)
        }
    }
    
    # one AAStringSet object:
    
    many_slices <- function(n) {
        sl <- unlist(sapply(
            X = unlist(smi[n]), FUN = slice, 
            seq = input_obj[n]
        ))
    }
    
    smt <- sapply(1:length(smi), many_slices)
    return(do.call(c, unlist(smt)))
    
}