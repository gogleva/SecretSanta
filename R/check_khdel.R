#' check C-terminal ER-retention signals
#'
#' This function checks the presence of terminal \strong{KDEL/HDEL} motifs in 
#' the provided amino acid equences.
#' @param input_obj input object of CBSResult class
#' @param pattern ER-pattern to check for:
#' \itemize{
#' \item \code{pattern = 'prosite'} - "[KRHQSA][DENQ]EL$>"
#' \item \code{pattern = 'elm'} - "[KRHQSAP][DENQT]EL$"  
#' \item \code{pattern = 'strict'} - "[KH]DEL$"
#' }
#' @return ErResult object
#' @export
#' @examples 
#' # check ER retention signals in CBSResult object
#' # before running signalp or any other
#' # predictions:
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta",
#'                                     package = "SecretSanta"))
#' inp <- CBSResult(in_fasta = aa[1:10])
#' # generate signalp predictions:
#' step1_sp2 <- signalp(inp, version = 4, organism = 'euk', 
#'     run_mode = 'starter')
#' # check ER retention signal in the signalp output, PROSITE pattern
#' et_result <- check_khdel(step1_sp2, pattern = 'prosite')

check_khdel <- function(input_obj, pattern = c('prosite', 'elm', 'strict')) {
    
    # check the input
    if (is(input_obj, "CBSResult")) {} else {stop(
        'input_object does not belong to CBSResult superclass')
    }
    
    if (length(getOutfasta(input_obj)) == 0) {
        stop('the input object contains empty out_fasta slot')
    }
    
    # check pattern argument
    
    # check organism argument
    if (missing(pattern)) {
        stop('missing argument: pattern')
    }
    pattern <- match.arg(pattern)
    
    # starting message:
    message(paste("checking for terminal ER retention signals,", pattern, 'pattern'))
    
    # select out_fasta slot
    fasta <- getOutfasta(input_obj)
    
    message(paste('Submitted sequences...', length(fasta)))
    
    # find and remove termial ER retention motifs
    
    prosite_er <- '[KRHQSA][DENQ]EL$'
    elm_er <- '[KRHQSAP][DENQT]EL$'
    strict_er <- '[KH]DEL$'
    
    # determine which pattern to check for
    
    if (pattern == 'prosite') {
        check_pattern <- prosite_er
    } else if (pattern == 'elm') {
        check_pattern <- elm_er
    } else if (pattern == 'strict') {
        check_pattern <- strict_er
    }    
    
    to_retain <- as.logical(lapply(fasta, function(x) grepl(check_pattern, x)))

    # generate cropped names for input fasta
    cropped_names <- unname(sapply(names(fasta), crop_names))
    # replace long names with cropped names
    names(fasta) <- cropped_names
    
    out_obj <- ErResult(retained_fasta = fasta[to_retain])
    out_obj <- setInfasta(out_obj, in_fasta = fasta)
    out_obj <- setOutfasta(out_obj, out_fasta = fasta[!to_retain])
    
    ret_count <- sum(to_retain)
    non_ret_count <- sum(!to_retain)
    message(paste('Sequences with terminal ER retention signals detected...',
                                ret_count))
    message(paste('Candidate without terminal ER retention signals detected...',
                                non_ret_count))
    
    if (validObject(out_obj)) {return(out_obj)}
}

