# HELPER FUNCTIONS FOR PARALLEL SIGNALP:

#' split XStringSet objects
#'
#' This function splits large XStringSet objects into chunks of given size and
#' returns a list of AAStringSet objects.
#' @param string_set input AAStringSet object;
#' @param chunk_size the number of sequenses in a single chunk;
#' @return list of AAStringSet chunks.
#' @export
#' @examples
#' # Read fasta file:
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta",
#' package = "SecretSanta"))
#' # Split it into chunks
#' # with 10 sequences each:
#' split_XStringSet(aa,10)

split_XStringSet <- function(string_set, chunk_size) {
    if (!(is(string_set, 'XStringSet'))) {
        stop('Input string_set does not belong to XStringSet class')
    }
    lss <- length(string_set)
    if (chunk_size > lss) {stop('Chunk size exceeds total seq number')}

    total_seq  <- c(1:lss)
    chunks <- split(total_seq,
                    ceiling(seq_along(total_seq) / chunk_size))

    lapply(chunks, function(x) string_set[x])
}

#' combine multiple objects of SignalpResult class
#'
#' This function combines multiple instances of SignalpResult class,
#' typically generated with parLapply while running signalp predictions in
#' parallel mode.
#' @param arguments a list of SignalpResult objects to be combined
#' @export
#' @return SignalpResult object
#' @examples
#' aa <- readAAStringSet(system.file(
#' "extdata", "sample_prot_100.fasta", package = "SecretSanta"))
#' inp2 <- CBSResult(in_fasta = aa[1:10])
#' inp3 <- CBSResult(in_fasta = aa[20:30])
#' inp4 <- CBSResult(in_fasta = aa[40:50])
#' r1 <- signalp(input_obj = inp2, version = 4, organism = 'euk',
#' run_mode = 'starter')
#' r2 <- signalp(input_obj = inp3, version = 4, organism = 'euk',
#' run_mode = 'starter')
#' r3 <- signalp(input_obj = inp4, version = 4, organism = 'euk',
#' run_mode = 'starter')
#' obj <- list(r1, r2, r3)
#' combined_sp <- combine_SpResult(obj)


combine_SpResult <- function(arguments) {
    if ((all(sapply(arguments, is, 'SignalpResult'))) == FALSE) {
    stop('Some objects from argument list do not belong to SignalpResult class')
    }

    if (length(unique(sapply(arguments, getSPversion))) != 1) {
    stop('Supplied objects were generated with different versions of signalp')
    }

    # combine accesors for fasta slots in a list
    fasta_getters <- list(getInfasta, getOutfasta, getMatfasta)

    # apply the list with fasta getters to all the elements in the argument list
    z <- Map(function(x) sapply(fasta_getters, function(f) f(x)), arguments)
    
    # combine fasta fields by slot types
    got_fastas <- sapply(1:3, 
                         function(i) do.call(c, sapply(z, function(x) x[i])))

    # arrange combined slots in SignalpResult object:
    c_obj <- SignalpResult(
        in_fasta = got_fastas[[1]],
        out_fasta = got_fastas[[2]],
        mature_fasta = got_fastas[[3]],
        sp_tibble = do.call(rbind, (lapply(arguments, getSPtibble))),
        sp_version = getSPversion(arguments[[1]])
    )
}

# PARALLEL SIGNALP ITSELF:

#' predict signal peptides
#'
#' This function calls the command line tool \code{signalp} to predict the presence and location of
#' signal peptide cleavage sites in amino acid sequences.
#' \cr
#' \cr
#' Large input files (>500 sequnces) are automatically split into smaller chunks
#' so that \code{signalp} prediction can be run as an embarassingly parallel process
#' on a specified number of cores.
#' @param input_obj   an instance of class \code{CBSResult} containing protein
#' sequences as one of the attributes
#' @param version version of \code{signalp}, supported versions include: \cr
#'  \code{signalp2}, \code{signalp3}, \code{signalp4}.
#' @param organism a character string with the following options:
#' \itemize{
#' \item \code{organism = "euk"} - for eukaryotes
#' \item \code{organism = "gram+"} - for gram-positive bacteria
#' \item \code{organism = "gram-"} - for gram-negative bacteria
#' }
#' @param run_mode a character string with the following options:
#' \itemize{
#' \item \code{run_mode = "starter"} - if it is the first step in pipeline
#' \item \code{run_mode = "piper"} - if you run this function on the output of other CBS tools
#' }
#' @param paths if required version of \code{signalp} is not acessible globally, a file
#' conatining a full path to it's executable should be provided; for details
#' please check SecretSanta vignette.
#' @param truncate a logical indicating:
#' \itemize{
#' \item \code{truncate = TRUE} - sequences longer 2000 residues will be
#' truncated to this length limit and renamed
#' \item \code{truncate = FALSE} - long sequences will be excluded from the analysis
#' }
#' Default is \code{truncate = TRUE}.
#' @param cores number of cores for multicore processing. Default is \code{cores = 1}.
#' @param legacy_method optional argument, which prediction method to use when running SiganlP 2.0 and SignalP 3.0:
#' \itemize{
#' \item \code{legacy_method = "hmm"} - for HMM-based predictions
#' \item \code{legacy_method = "nn"} - for prediction based on neural networks
#' }
#' @return an object of SignalpResult class
#' @export
#' @examples
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta",
#' package = "SecretSanta"))
#' # assign this object to the input_fasta slot
#' # of empty CBSResult object
#' inp <- CBSResult(in_fasta = aa[1:10])
#' # run signalp2 on the initial file:
#' r1 <- signalp(inp, version = 2, organism = 'euk', run_mode = "starter",
#' legacy_method = 'hmm')
#' r4 <- signalp(inp, version = 4, organism = 'euk', run_mode = "starter")
#' @seealso \code{\link{parse_signalp}}

signalp <- function(input_obj,
                    version,
                    organism = c('euk', 'gram+', 'gram-'),
                    run_mode = c('starter', 'piper'),
                    paths = NULL,
                    truncate = NULL,
                    cores = 1,
                    legacy_method) {

  if (!is.numeric(cores))
    stop('cores argument must be numeric', call. = FALSE)

  if (parallel::detectCores() < cores)
    stop(
      "your machine has only ",
      parallel::detectCores() ,
      " cores ... you specified: ",
      cores,
      " cores. please fix ...",
      call. = FALSE
    )
    
    # ----- Check that inputs are valid
    # arguments are present and have valid values:
    if (missing(organism)) {
        stop('missing argument: organism', call. = FALSE)
    }

    if (!is.element(organism, c("euk", "gram+", "gram-")))
      stop("please specify an organism that is supported by this function ...", call. = FALSE)

    if (missing(run_mode)) {
        stop('missing argument: run_mode', call. = FALSE)
    }
    
    if (!is.element(run_mode, c("starter", "piper")))
      stop("please specify a run_mode that is supported by this function ...", call. = FALSE)
    organism <- match.arg(organism)
    run_mode <- match.arg(run_mode)

    
    # check that input object belongs to CBSResult class
    if (is(input_obj, "CBSResult")) {
    } else {
        stop('input_object does not belong to CBSResult superclass')
    }

    # check that input_object contains non-empty in/out_fasta for starter/piper
    if (run_mode == 'starter') {
        if (length(getInfasta(input_obj)) != 0) {
            fasta <- getInfasta(input_obj)
        } else {
            stop('in_fasta attribute is empty ...', call. = FALSE)
        }

    } else if (run_mode == 'piper') {
        if (length(getOutfasta(input_obj)) != 0) {
            fasta <- getOutfasta(input_obj)
        } else {
            stop('out_fasta attribute is empty ...', call. = FALSE)
        }
    }

    # check that version number is valid:
    if (is.element(version, c(2, 3, 4))) {
    } else {
        stop('version is invalid, allowed versions: c(2,3,4)', call. = FALSE)
    }

    # ----- Set default value for parameters if not provided:

    if (is.null(truncate))
        truncate <- TRUE
    else
        truncate
    if (!is.logical(truncate))
      stop("truncate can only be TRUE or FALSE ...", call. = FALSE)

    # ------ Produce encouraging status messages, inputs should be ok.
    # create actual tool name with version number provided
    signalp_version <- paste("signalp", version, sep = '')
    message(paste('Version used ...', signalp_version))
    message("running signalp locally ...")

    # simple signalp, takes single AAStringSet as an input and runs
    # signalp prediction on it

    simple_signalp <- function(aaSet) {
        # ---- Run prediction

        # convert fasta to a temporary file:
        out_tmp <- tempfile() #create a temporary file for fasta
        writeXStringSet(aaSet, out_tmp) #write tmp fasta file to use later
        message(paste('Submitted sequences...', length(aaSet)))

        # get and check paths to signalp
        if (is.null(paths)) {
            full_pa <- signalp_version
        } else {
            mp <- suppressMessages(manage_paths(
                in_path = FALSE,
                test_tool = signalp_version,
                path_file = paths))
            full_pa <- mp$path_tibble$path
        }

        # decide how to run signalp depending on a version provided

        if (version == 4) {
            # runing signalp versios 4 and 4.1, potentially should work for 5

            sp <-
                tibble::as.tibble(read.table(text = (system(
                    paste(full_pa, "-t", organism,
                          "-f short -trunc 70", out_tmp),
                    intern = TRUE
                ))))
            names(sp) <- c("gene_id", "Cmax", "Cpos",
                           "Ymax", "Ypos", "Smax",
                           "Spos", "Smean", "D",
                           "Prediction", "Dmaxcut",
                           "Networks-used")

            # reorder columns to match sp2/3 output:
            sp <- sp[c("gene_id", "Cmax", "Cpos",
                       "Ymax", "Ypos", "Smax",
                       "Spos", "Smean", "Prediction")]

            sp <- sp[sp$Prediction == 'Y',]
            sp$Prediction <- ifelse(sp$Prediction == 'Y', 'Signal peptide', "none")

        } else if (version < 4) {
            # running signalp versions 2 and 3, call parse_signalp
            message('signalp < 4, calling parser for the output...')
            con <-
                system(paste(full_pa, "-t", organism, "-f short -trunc 70",
                             "-m", legacy_method, out_tmp), intern = TRUE)

            sp <- parse_signalp(input = con,
                                input_type = "system_call",
                                version = version,
                                method = legacy_method,
                                source_fasta = out_tmp,
                                pred_filter = "Signal peptide")
        }

        message(paste('Candidate sequences with signal peptides...', nrow(sp)))

        if (nrow(sp) == 0) {
            warning('Signal peptide prediction yeilded 0 candidates ... ', call. = FALSE)
        }

        # generate cropped names for input fasta
        names(aaSet) <- unname(sapply(names(aaSet), crop_names))
        # get ids of candidate secreted proteins
        out_fasta_sp <- aaSet[sp$gene_id]

        # generate mature sequences
        cropped_fasta <- subseq(out_fasta_sp, start = sp$Cpos, end = -1)

        # construct output object
        out_obj <- SignalpResult(
            in_fasta = aaSet,
            out_fasta = out_fasta_sp,
            mature_fasta = cropped_fasta,
            sp_version = version,
            sp_tibble = sp
        )

        # check that intended output is valid
        if (validObject(out_obj)) {
            return(out_obj)
        }
    } #close simple_fasta

    # Handle long sequences if any present in the input,
    # if necessary - run signalp as a parallel process

    # estimate how big is the file, if required - split it into smaller
    # chunks and run signalp as an embarassingly parallel process

    fasta <- truncate_seq(truncate = truncate, fasta, 2000)
    if (length(fasta) <= 600) {
        message('Ok for single processing')
        return(simple_signalp(estimate_lim(fasta, truncate)))
    } else {
        message('Input fasta contains >600 sequences, entering batch mode...')

        #split and check that chunks do not exceed 200K residue limit
        split_fasta <-
            sapply(split_XStringSet(fasta, 500), estimate_lim, truncate)

        # Initiate cluster
        cl <- parallel::makeCluster(cores)
        # run parallel process

        parallel::clusterExport(cl = cl, varlist = c("paths"), envir = environment())
        result <- parallel::parLapply(cl, split_fasta, simple_signalp)
        parallel::stopCluster(cl)

        res_comb <- do.call(c, result)
        combined_SignalpResult <- combine_SpResult(unname(res_comb))

        sp_count <- nrow(getSPtibble(combined_SignalpResult))
        message(paste('Candidate sequences with signal peptides...', sp_count))
        if (sp_count == 0) {
            warning('Signal peptide prediction yeilded 0 candidates ...', call. = FALSE)
        }

    closeAllConnections()
    return(combined_SignalpResult)
    }
}
