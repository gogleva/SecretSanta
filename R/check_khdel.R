#' check C-terminal ER-retention signals
#'
#' This function checks the presence of terminal \strong{KDEL/HDEL} motifs in 
#' the provided amino acid equences.
#' @param input_obj input object of CBSResult class
#' @return ErResult object
#' @export
#' @examples 
#' # check ER retention signals in CBSResult object
#' # before running signalp or any other
#' # predictions:
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta",
#'                                     package = "SecretSanta"))
#' inp <- SignalpResult(in_fasta = aa[1:10])
#' et_s <- check_khdel(inp,
#'     run_mode = 'starter')
#' 
#' # generate signalp predictions:
#' step1_sp2 <- signalp(et_s, version = 4, organism = 'euk', 
#'     run_mode = 'starter')
#' # check ER retention signal in the signalp output
#' et_result <- check_khdel(step1_sp2)

check_khdel <- function(input_obj) {
    
    # check the input
    if (is(input_obj, "CBSResult")) {} else {stop(
        'input_object does not belong to CBSResult superclass')
    }
    
    if (length(getOutfasta(input_obj)) == 0) {
        stop('the input object contains empty out_fasta slot')
    }
    
    # starting message:
    message("checking for terminal ER retention signals...")
    
    # select out_fasta slot
    fasta <- getOutfasta(input_obj)
    
    message(paste('Submitted sequences...', length(fasta)))
    
    # find and remove termial ER retention motifs
    ER1 <- Biostrings::AAString('KDEL')
    ER2 <- Biostrings::AAString('HDEL')
    tails <- subseq(fasta, -4, -1) # last 4 aminno acids in each protein
    un <- !(as.logical(vcountPattern(ER1, tails)) |
                        as.logical(vcountPattern(ER2, tails)))
    
    # generate cropped names for input fasta
    cropped_names <- unname(sapply(names(fasta), crop_names))
    # replace long names with cropped names
    names(fasta) <- cropped_names
    
    # fasta without terminal KDELs/HDELs
    non_retained <- fasta[un] 
    
    out_obj <- ErResult(in_fasta = fasta,
                        out_fasta = non_retained,
                        retained_fasta = fasta[!un])
    
    ret_count <- length(fasta[!un])
    message(paste('Sequences with terminal ER retention signals detected...',
                                ret_count))
    message(paste('Candidate without terminal ER retention signals detected...',
                                length(non_retained)))
    
    if (validObject(out_obj)) {return(out_obj)}
}
