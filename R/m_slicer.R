#' generate proteins with alternative translation start sites
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
#' @param run_mode    \strong{slice} - just slice input fasta, regardless
#' of its origin; \cr
#' \strong{rescue} - get proteins not predicted to be secreted on the 
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

#FOR_TESTING----
#input_obj <- readAAStringSet(system.file("extdata","sample_prot_100.fasta",
#                                        package = "SecretSanta"))
#input_obj <- CBSResult(in_fasta = aa[1:10])
###FOR_TESTING----

m_slicer <- function(input_obj, min_len, run_mode = c('slice', 'rescue')) {

    # check that inputs are present and valid
    if (missing(run_mode)) {stop('missing argument: run_mode')}
   
    run_mode = match.arg(run_mode)
    
    if (is(input_obj, 'AAStringSet') && (run_mode != 'slice')) {
        stop("Use run_mode 'slice' for an input object of AAStringSet class")
        }
    
    if (is(input_obj, 'CBSResult') && (run_mode != 'rescue')) {
    stop("Use run_mode 'rescue' for an input object of CBSResult class")
        }
    
    # transform input_object if necessary
    
    if (is(input_obj, 'CBSResult')) {
        infa <- getInfasta(input_obj)
        outfa <- getOutfasta(input_obj)
        input_obj <- infa[!is.element(names(infa), names(outfa))]
    } else {input_obj}
    
    # get positions of all Meths in the AAStringSet object
    
    mi <- vmatchPattern('M', input_obj)
    
    # select start M-positions and remove first methionine from all M-coord 
    # lists (1st, as in the original protein sequences)
    smi <- sapply(startIndex(mi), function(x) x[x > 1])
    
    # filter input AAStringSet object and M start indexes list (to remove 
    # sequences that don't have alteranative methionines, probbly will be none,
    # but it is better to check
    
    input_obj <- input_obj[lengths(smi) > 0]
    smi <- smi[lengths(smi) > 0]
    
    # function to slice one AAString based on single M coordinate
    slice_string <- function(x, seq) {
        st <- subseq(seq, start = x, end = -1)
        # cop and update names by adding with M-coordinates sequences were
        # sliced from
        names(st) <- paste(strsplit(names(st), ' ')[[1]][1],
                        '_slice_M',
                        x,
                        sep = '')
        if (width(st) >= min_len) return(st)
        }
    
    # extention to slice AAStringSet (multiple strings):
    # may be just an anonymous function instead?
    # do I even need this unlist?
    
    slice_set <- function(n) {
        unlist(sapply(smi[[n]], slice_string, input_obj[n]))
    }
    
    smt <- sapply(1:length(smi), slice_set)
    return(do.call(c, unlist(smt)))
    
}
