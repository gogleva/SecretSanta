#' predict transmembrane domains with TMHMM
#'
#' This function calls local TMHMM to predict transmembarne domains in protein
#' sequence. \cr
#' \cr
#' In the full-length proteins N-terminal signal peptide could be
#' erroneously predicted as a TM domain. Therefore it is recommended to run
#' TMHMM on mature pepetide sequences (with clipped signal peptides) to avoid 
#' false positive predictions. \cr
#' \cr
#' To generate an CBSResult object with mature sequences please run signalp
#' predictions first.
#' @param input_obj an instance of SignalpResult class, \cr
#'                  input should contain mature_fasta slot;
#' @param TM  allowed number of TM domains in mature peptides,
#' recommended value <= 1; use TM = 0 for strict filtering                   
#' @param paths if tmhmm is not acessible globally, a file
#' conatining a full path to it's executable should be provided; for details
#' please check SecretSanta vignette.
#' @export
#' @return TMhmmResult object
#' @examples 
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta",
#' package = "SecretSanta"))
#' inp <- CBSResult(in_fasta = aa[1:10])
#' s1_sp2 <- signalp(inp, version = 2, organism = 'euk', run_mode = "starter")
#' tm <- tmhmm(s1_sp2, TM = 1)

tmhmm <- function(input_obj, TM, paths = NULL) {
    
    # ----- Check the inputs
    
    if (is.numeric(TM)) {
    } else {
        stop('TM argument should be numeric')
    }
    
    if (TM >= 2) {
        warning('Recommended TM threshold values for mature peprides is 1')
    }
    
    # check that input object belongs to a valid class
    if (is(input_obj, "SignalpResult")) {
    } else {
        stop('input_object does not belong to SignalpResult class')
    }
    
    # check that input object contains non-empty mature fasta slot
    s <- getSlots(class(input_obj))
    
    if ('mature_fasta' %in% names(s)) {
        if (length(getMatfasta(input_obj)) == 0) {
            stop('the input object contains an empty mature_fasta slot')
        }
    } else {
        stop('the input object does not contain mature_fasta slot')
    }
    
    #----- Run tmhmm
    message("running TMHMM locally...")
    
    fasta <- getMatfasta(input_obj)
    out_tmp <- tempfile()
    Biostrings::writeXStringSet(fasta, out_tmp)
    
    message(paste('Submitted sequences...', length(fasta)))
    
    # get and check paths to tmhmm
    if (is.null(paths)) {
        full_pa <- 'tmhmm'
    } else {
        mp <- suppressMessages(manage_paths(
            in_path = FALSE,
            test_tool = 'tmhmm',
            path_file = paths
        ))
        full_pa <- mp$path_tibble$path
    } 

    con <- system(paste(full_pa, out_tmp, '--short'), intern = TRUE)
    con_tmp <- tempfile()
    write(con, con_tmp)
    tm <- tibble::as.tibble(read.table(con_tmp, sep = '\t', header = FALSE,
                                       stringsAsFactors = FALSE))
    
    names(tm) <- c("gene_id", "length", "ExpAA", "First60", "PredHel",
                    "Topology")
    
    # clean output values remove '... =' value
    clean_outp <- function(x) {
        unlist(strsplit(x, '='))[2]
    }
    tm <- dplyr::mutate(
        tm,
        length = as.numeric(sapply(tm$length, clean_outp, USE.NAMES = FALSE)),
        ExpAA = as.numeric(sapply(tm$ExpAA, clean_outp, USE.NAMES = FALSE)),
        First60 = as.numeric(sapply(tm$First60, clean_outp, USE.NAMES = FALSE)),
        PredHel = as.numeric(sapply(tm$PredHel, clean_outp, USE.NAMES = FALSE)),
        Topology = sapply(tm$Topology, clean_outp, USE.NAMES = FALSE)
    )
    
    # filter based on TM threshold
    tm <- (tm %>% dplyr::filter_( ~ PredHel <= TM))
    
    message(paste(
        'Candidates with signal peptides and 0 TM domains in mature seq',
        nrow(tm)
    ))
    
    #generate cropped names for input fasta
    full_fasta <- getInfasta(input_obj)
    cropped_names <- unname(sapply(names(full_fasta), crop_names))
    
    #replace long names with cropped names
    names(full_fasta) <- cropped_names
    
    #get ids of candidate secreted proteins
    candidate_ids <- tm %>%
        dplyr::select_( ~ gene_id) %>%
        unlist(use.names = FALSE)
    out_fasta_tm <- full_fasta[candidate_ids]
    
    out_obj <- TMhmmResult(
            in_mature_fasta = fasta,
            out_mature_fasta = fasta[candidate_ids],
            tm_tibble = tm)
    # original in fasta
    out_obj <- setInfasta(out_obj, in_fasta = getOutfasta(input_obj))
    # out fasta, full length
    out_obj <- setOutfasta(out_obj,  out_fasta = out_fasta_tm)
    # clean TMP files before exiting:
    
    junk <- dir(pattern = 'TMHMM*')
    file.remove(junk)
    
    if (validObject(out_obj)) {
        return(out_obj)
    }
}
