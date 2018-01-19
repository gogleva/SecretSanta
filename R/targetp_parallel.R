#' combine multiple objects of TargetpResult class
#'
#' This function combines multiple instances of TargetpResult class,
#' typically generated with parLapply when running targetp in parallel mode.
#' @param arguments  a list of TargetpResult objects to be combined.
#' @export
#' @return TargetpResult object
#' @examples 
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta",
#' package = "SecretSanta"))
#' inp2 <- CBSResult(in_fasta = aa[1:10])
#' inp3 <- CBSResult(in_fasta = aa[20:30])
#' inp4 <- CBSResult(in_fasta = aa[40:50])     
#' tp1 <- targetp(input_obj = inp2, network = 'N', run_mode = 'starter')
#' tp2 <- targetp(input_obj = inp3, network = 'N', run_mode = 'starter')
#' tp3 <- targetp(input_obj = inp4, network = 'N', run_mode = 'starter')
#' obj <- list(tp1, tp2, tp3)
#' combined_tp <- combine_TpResult(obj)

combine_TpResult <- function(arguments) {
    if (all(sapply(arguments, is, 'TargetpResult'))) {
        
    } else {
    stop('Some objects from argument list do not belong to TargetpResult class.')
    }
    
    c_in_fasta <- do.call(c, (lapply(arguments, getInfasta)))
    c_out_fasta <- do.call(c, (lapply(arguments, getOutfasta)))
    c_tp_tibble <- do.call(rbind, (lapply(arguments, getTPtibble)))
    
    c_obj <- TargetpResult(tp_tibble = c_tp_tibble)
    c_obj <- setInfasta(c_obj, in_fasta = c_in_fasta)
    c_obj <- setOutfasta(c_obj, out_fasta = c_out_fasta)
    
}

#' predict subcellular protein localization with TargetP
#'
#' This function calls local targetp to predict subcellular localisation of 
#' a protein.\cr
#' \cr
#' For large inputs (\strong{>1000} sequences) will automatically run
#' as a massive paralle job.
#' @param input_obj    an instance of CBSResult class containing protein 
#' sequences in one of the slots
#' @param network    
#' \strong{P} - for plants; \cr
#' \strong{N} - for non-plants;
#' @param run_mode
#' \strong{starter} - if it is the first step in pipeline; \cr
#' \strong{piper} - if you run this function on the output of other CBS tools;
#' @param paths if targetp is not acessible globally, a file
#' conatining a full path to it's executable should be provided; for details
#' please check SecretSanta vignette.
#' @param cores optional arguments, number of cores to run the parallel process
#' on. If not set default will be 1.  
#' @export
#' @return an object of TargetpResult class
#' @examples 
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta",
#' package = "SecretSanta"))
#' # assign this object to the input_fasta 
#' # slot of CBSResult object
#' inp <- CBSResult(in_fasta = aa[1:10])
#' # run target prediction
#' tp_result <- targetp(input_obj = inp, network = 'N', run_mode = 'starter')

targetp <- function(input_obj,
                    network = c('P', 'N'),
                    run_mode = c('starter', 'piper'),
                    paths = NULL,
                    cores = NULL) {
    
    # ----- Check that inputs are valid
    
    # check that arguments are present and valid
    if (missing(network)) {
        stop('Missing argument: network.')
    }
    if (missing(run_mode)) {
        stop('Missing argument: run_mode.')
    }
    organism <- match.arg(network)
    run_mode <- match.arg(run_mode)
    
    if (is.null(cores))
        cores = 1
    else
        cores
    if (is.numeric(cores)) {
    } else {
        stop('Cores argument must be numeric.')
    }
    if (cores > detectCores()) {
        stop('Cores value > available core number.')
    }
    
    # check that input object belong to CBSResult class
    if (is(input_obj, "CBSResult")) {
    } else {
        stop('input_obj does not belong to CBSResult superclass.')
    }
    
    # check that input_obj contains non-empty in/out_fasta for starter/piper
    if (run_mode == 'starter') {
        if (length(getInfasta(input_obj)) != 0) {
            fasta <- getInfasta(input_obj)
        } else {
            stop('in_fasta attribute is empty.')
        }
    } else if (run_mode == 'piper') {
        if (length(getOutfasta(input_obj)) != 0) {
            fasta <- getOutfasta(input_obj)
        } else {
            stop('out_fasta attribute is empty.')
        }
    }
    
    # All checked, produce an encouragig message
    message("running targetp locally...")
    
    # simple function to run targetp on small input files <1K proteins
    
    simple_targetp <- function(aaSet) {
        #----- Run targetp prediction:
        message(paste('Number of submitted sequences...', length(aaSet)))
        
        # convert fasta to a temporary file:
        out_tmp <- tempfile() #create a temporary file for fasta
        Biostrings::writeXStringSet(aaSet, out_tmp) #write tmp fasta file
        
        # get and check paths to signalp
        if (is.null(paths)) {
            full_pa <- 'targetp'
        } else {
            mp <- suppressMessages(manage_paths(
                in_path = FALSE,
                test_tool = 'targetp',
                path_file = paths
            ))
            full_pa <- mp$path_tibble$path
        }
        
        # prep fasta:
        # generate cropped names for input fasta
        cropped_names <- unname(sapply(names(aaSet), crop_names))
        # replace long names with cropped names
        names(aaSet) <- cropped_names
        
        # prep networks argument:
        NN <- paste('-', network, sep = '')
        
        #run targetp:
        tp <- tibble::as.tibble(read.table(text = (system(
            paste(full_pa, NN, out_tmp),
            intern = TRUE
        )[1:length(aaSet) + 8])))
        
        if (network == 'N') {
            names(tp) <- c('gene_id', 'length', 'mTP', 'sp', 'other',
                            'TP_localization', 'RC')
        }
        else if (network == 'P') {
            names(tp) <- c('gene_id', 'length', 'cTP', 'mTP', 'sp', 'other',
                            'TP_localization', 'RC')
        }
        
        # replace gene_id with the full ones, TargetP crops everything to 20 characters
        tp$gene_id <- names(aaSet)
        tp <- tp %>% dplyr::filter_( ~ TP_localization == 'S')
        message(paste('Number of candidate secreted sequences', nrow(tp)))
        
        candidate_ids <- tp %>%
            dplyr::select_( ~ gene_id) %>%
            unlist(use.names = FALSE)
        out_fasta_tp <- aaSet[candidate_ids]
        
        # generate output object:
        out_obj <- TargetpResult(tp_tibble = tp)
        out_obj <- setInfasta(out_obj, in_fasta = aaSet)
        out_obj <- setOutfasta(out_obj, out_fasta = out_fasta_tp)
        
        if (validObject(out_obj)) {
            return(out_obj)
        }
    }

    # Check input file size and decide how to run targetp:
    # in parallel mode or not:
    
    if (length(fasta) <= 1000) {
        message('Ok for single processing')
        return(simple_targetp(fasta))
    } else {
        message('Input fasta contains >1000 sequences, entering batch mode...')
        message(paste('Number of submitted sequences...', length(fasta)))
        # split fasta:
        split_fasta <- split_XStringSet(fasta, 1000)
        
        # initiate cluster
        cl <- makeCluster(cores)
        # clusterEvalQ(cl, library(SecretSanta))
        
        clusterExport(cl = cl, varlist = c("paths"), envir = environment())
        
        # run parallel targetp:
        result <- parLapply(cl, split_fasta, simple_targetp)
        stopCluster(cl)
        
        # combine outputs from multiple workers
        res_comb <- do.call(c, result)
        combined_TargetpResult <- combine_TpResult(unname(res_comb))
        
        tp_count <- nrow(getTPtibble(combined_TargetpResult))
        message(paste('Number of candidate secerted sequences...', tp_count))
        if (tp_count == 0) {
            warning('Targetp prediction yeilded 0 extracellular candidates')
        }
        
        closeAllConnections()
        return(combined_TargetpResult)
        
    }
}


