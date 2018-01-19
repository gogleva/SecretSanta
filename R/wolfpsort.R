#' predict subcellular protein localization with WoLFpsort
#'
#' This function runs WoLF PSORT to predict protein cellular sub-localisation
#' and returns the most probbale one. Including this step in secretome 
#' prediction pipelines provides additional supportig evidence that a protein
#' might be secreted and deposited outside the cell.\cr
#' \cr
#' Recommended to run on the late stages of secretome prediction pipeline.\cr
#' \cr
#' Also see targetp function - for similar functionality.
#' 
#' @param input_obj Object of CSBResult class
#' @param organism  set relevant taxonomic group,
#'                  options include: \strong{plant},
#'                  \strong{animal}, \strong{fungi};
#' @param paths   if wolfpsort is not acessible globally, a file
#' conatining a full path to it's executable should be provided; for details
#' please check SecretSanta vignette. 
#' @param run_mode 
#' \strong{starter} - if it is the first step in pipeline; \cr
#' \strong{piper} - if you run this function on the output of other methods;
#' @return object of WolfResult class  
#' @export
#' @examples
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata","sample_prot_100.fasta",
#' package = "SecretSanta"))
#' # assign this object to the input_fasta slot of CBSResult object
#' inp <- CBSResult(in_fasta = aa[1:10])
#' # run signalp2 on the initial file:
#' step1_sp2 <- signalp(inp, version = 2, organism ='euk',
#' run_mode = "starter", legacy_method = 'hmm')
#' # run wolfpsort on the signalp output:
#' w <- wolfpsort(step1_sp2, 'fungi', run_mode = 'piper')

wolfpsort <- function(input_obj, organism = c('plant', 'animal', 'fungi'),
                      run_mode = c('starter', 'piper'), paths = NULL) {
    # check that inputs are valid
    
    # check organism argument
    if (missing(organism)) {
        stop('Missing argument: organism.')
    }
    organism <- match.arg(organism)
    
    # check input_obj
    if (is(input_obj, "CBSResult")) {
    } else {
        stop('input_object does not belong to CBSResult superclass.')
    }

    message("Running WoLF PSORT locally ...")
    
   # fasta <- getOutfasta(input_obj)

   # check that input_obj contains non-empty in/out_fasta for starter/piper
    if (run_mode == 'starter') {
        if (length(getInfasta(input_obj)) != 0) {
            fasta <- getInfasta(input_obj)
        } else {
            stop('in_fasta slot is empty.')
        }
    } else if (run_mode == 'piper') {
        if (length(getOutfasta(input_obj)) != 0) {
            fasta <- getOutfasta(input_obj)
        } else {
            stop('out_fasta slot is empty')
        }
    }
   
    message(paste("Number of submitted sequences ...", length(fasta)))
    
    # prep fasta:
    # generate cropped names for input fasta
    cropped_names <- unname(sapply(names(fasta), crop_names))
    # replace long names with cropped names
    names(fasta) <- cropped_names
    
    out_tmp <- tempfile()
    Biostrings::writeXStringSet(fasta, out_tmp)
    
    # get and check paths to wolfpsort
    if (is.null(paths)) {
        full_pa <- 'wolfpsort'
    } else {
        mp <- suppressMessages(manage_paths(
            in_path = FALSE,
            test_tool = 'wolfpsort',
            path_file = paths
        ))
        full_pa <- mp$path_tibble$path
    }

    wolf <-
        system(paste(full_pa,  organism, '<', out_tmp), intern = TRUE)
    
    #parse wolf output
    clean_strings <- function(x, field) {
        unlist((strsplit(x, " ")))[c(field)]
    }
    
    gene_id <- sapply(X = wolf, field = 1, clean_strings, USE.NAMES = FALSE)
    localization <- sapply(X = wolf, field = 2, clean_strings, 
                            USE.NAMES = FALSE)
    
    #assemble result tibble with gene id and most probable
    #subsellular localisation
    wolf_tbl <-
        tibble::as_tibble(data.frame(gene_id, localization)) %>%
        filter_( ~ gene_id != '#') %>%
        filter_( ~ localization == 'extr')
    
    message(paste(
        'Candidate sequences with extracellular localisation ...',
        nrow(wolf_tbl)
    ))
    
    #assemble wolf result object:
    out_obj <- WolfResult(wolf_tibble = wolf_tbl)
    out_obj <- setInfasta(out_obj, in_fasta = fasta)
    out_obj <- setOutfasta(out_obj, out_fasta = fasta[wolf_tbl$gene_id])
    
    if (validObject(out_obj)) {
        return(out_obj)
    }
}
