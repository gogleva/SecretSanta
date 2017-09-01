# HELPER FUNCTIONS FOR PARALLEL SIGNALP:

#' split_XStringSet function
#'
#' This function splits large XStringSet objects into chunks of given size and 
#' returns a list of AAStringSet objects.
#' @param string_set input AAStringSet object;
#' @param chunk_size the number of sequenses in a single chunk;
#' @return list of AAStringSet chunks.
#' @export
#' @examples 
#' # Read fasta file:
#' 
#' aa <- readAAStringSet(
#'  system.file(
#'    "extdata", 
#'    "sample_prot_100.fasta",
#'     package = "SecretSanta"
#'  ))
#'  
#' # Split it into chunks
#' # with 10 sequences each: 
#'                                    
#' split_XStringSet(aa,10)                                   

split_XStringSet <- function(string_set, chunk_size) {
  if (!(is(string_set, 'XStringSet'))) {
    stop('Input string_set does not belong to XStringSet class')
  }
  
  lst <- length(string_set)
  
  if (chunk_size > lst) {
    stop('Chunk size exceeds total seq number')
  }
  
  total_seq  <- c(1:lst)
  chunks <- split(total_seq,
                  ceiling(seq_along(total_seq) / chunk_size))
  
  seq_chunker <- function(x) {
    chunk <- string_set[x]
  }
  
  lapply(chunks, seq_chunker)
}

#' combine_SpResult function
#'
#' This function combines multiple instances of SignalpResult class,
#' typically generated with parLapply while running signalp predictions in 
#' parallel mode.
#' @param arguments a list of SignalpResult objects to be combined
#' @export
#' @return SignalpResult object
#' @examples 
#' aa <- readAAStringSet(system.file("extdata",
#'                  "sample_prot_100.fasta",
#'                   package = "SecretSanta"))
#' inp2 <- SignalpResult(in_fasta = aa[1:10])
#' inp3 <- SignalpResult(in_fasta = aa[20:30])
#' inp4 <- SignalpResult(in_fasta = aa[40:50])                   
#' r1 <- signalp(input_obj = inp2,
#'                version = 4,
#'                organism = 'euk',
#'                run_mode = 'starter')
#' r2 <- signalp(input_obj = inp3,
#'                version = 4,
#'                organism = 'euk',
#'                run_mode = 'starter')
#' r3 <- signalp(input_obj = inp4,
#'                version = 4,
#'                organism = 'euk',
#'                run_mode = 'starter') 
#' obj <- list(r1, r2, r3)
#' combined_sp <- combine_SpResult(obj)


combine_SpResult <- function(arguments) {
  if ((all(sapply(arguments, is, 'SignalpResult'))) == FALSE) {
    stop('Some objects from arguments list do not belong to SignalpResult class.')
  }
  
  c_in_fasta <- do.call(c, (lapply(arguments, getInfasta)))
  c_out_fasta <- do.call(c, (lapply(arguments, getOutfasta)))
  c_mature_fasta <- do.call(c, (lapply(arguments, getMatfasta)))
  c_sp_tibble <- do.call(rbind, (lapply(arguments, getSPtibble)))
  c_sp_version <- unlist((lapply(arguments, getSPversion))[1])
  
  c_obj <- SignalpResult(
    in_fasta = c_in_fasta,
    out_fasta = c_out_fasta,
    mature_fasta = c_mature_fasta,
    sp_tibble = c_sp_tibble,
    sp_version = c_sp_version
  )
}

# PARALLEL SIGNALP ITSELF:

#' signalp function
#'
#' This function calls local signalp to predict the presence and location of
#' signal peptide cleavage sites in amino acid sequences.
#' \cr
#' \cr
#' Large input files (>500 sequnces) are automatically split into smaller chunks
#' so that signalp prediction could be run as an embarassingly parallel process
#' on all the CPUs available.
#' @param input_obj   an instance of CBSResult class containing protein 
#' sequences as one of the attributes
#' @param version  signalp version to run, supported versions include: \cr
#'  2, 3, 4.
#' @param organism 
#' \strong{euk} - for eukaryotes;\cr
#' \strong{gram+} - for gram-positive bacteria;\cr
#' \strong{gram-} - for gram-negative bacteria;\cr
#' @param run_mode
#' \strong{starter} - if it is the first step in pipeline;\cr
#' \strong{piper} - if you run this function on the output of other CBS tools;
#' @param paths if required version of signalp is not acessible globally, a file
#' conatining a full path to it's executable should be provided; for details
#' please check SecretSanta vignette.
#' @param truncate if \strong{TRUE} - sequences longer 2000 residues will be 
#' truncated to this length limit and renamed;\cr
#' if \strong{FALSE} - long sequences will be excluded from the analysis;\cr
#' Default = TRUE.
#' @param cores optional arguments, number of cores to run the parallel process
#' on. If not set default will be 1.
#' @return an object of SignalpResult class
#' @export
#' @examples
#' 
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file(
#'                        "extdata",
#'                        "sample_prot_100.fasta",
#'                         package = "SecretSanta"))
#'                                   
#' # assign this object to the input_fasta slot
#' # of empty CBSResult object
#' inp <- CBSResult(in_fasta = aa[1:10])
#' 
#' # run signalp2 on the initial file:
#' r1 <- signalp(inp,
#'               version = 2,
#'               organism = 'euk'
#'               run_mode = "starter")

signalp <- function(input_obj,
                    version,
                    organism = c('euk', 'gram+', 'gram-'),
                    run_mode = c('starter', 'piper'),
                    paths = NULL,
                    truncate = NULL,
                    cores = NULL) {
  
  # ----- Check that inputs are valid
  
  # arguments are present and have valid values:
  if (missing(organism)) {stop('missing argument: organism')}
  if (missing(run_mode)) {stop('missing argument: run_mode')}
  organism <- match.arg(organism)
  run_mode <- match.arg(run_mode)
  
  # check that input object belong to CBSResult class
  if (is(input_obj, "CBSResult")) {} else {
    stop('input_object does not belong to CBSResult superclass')
  }
  
  # check that input_object contains non-empty in/out_fasta for starter/piper
  if (run_mode == 'starter') {
    if (length(getInfasta(input_obj)) != 0) {
      fasta <- getInfasta(input_obj)
    } else {stop('in_fasta attribute is empty')}
    
  } else if (run_mode == 'piper') {
    if (length(getOutfasta(input_obj)) != 0) {
      fasta <- getOutfasta(input_obj)
    } else {
      stop('out_fasta attribute is empty')}
    }
  
  # check that version number is valid:
  if (version %in% c(2,3,4)) {} else {
    stop('version is invalid, allowed versions: c(2,3,4)')}
  
  # ----- Set default value for parameters if not provided:
  
  if (is.null(truncate)) truncate = TRUE else truncate
  if (is.logical(truncate)) {} else {stop('truncate argument must be logical')}
  
  if (is.null(cores)) cores = 1 else cores
  if (is.numeric(cores)) {} else {stop('cores argument must be numeric')}
  if (cores > detectCores()) {stop('cores value > available core number')}

  # ------ Helper functions:
  
  # helper function: crop long names for AAStringSet object,
  # return character vector
  crop_names <- function(x) {
    unlist(stringr::str_split(x, " "))[1]
  }
  
  # helper function to truncate long sequences or throw them away, otherwise
  # signalp will break (at least signalp2 and signalp3 will)
  truncate_seq <- function(truncate, seq_set, threshold) {
    drop_n <- length(seq_set[width(seq_set) >= threshold])
    
    if (drop_n == 0) return(seq_set)
    
    if (truncate == FALSE) {
      seq_set <- seq_set[width(seq_set) < threshold]
      warning(paste(drop_n, 'long sequenses have been thrown away'))
      return(seq_set)
      
    } else if (truncate == TRUE) {
      message(paste(drop_n, 'sequences to be truncated'))
      seq_keep <- seq_set[width(seq_set) < threshold] # not so long sequences
      seq_trunc <- seq_set[width(seq_set) >= threshold] # sequences to truncate
      t_names <- paste(unname(sapply(names(seq_trunc), crop_names)),
                       '_truncated',
                       sep = '')
      names(seq_trunc) <- t_names #new names for sequences to be truncated
      seq_trunc <- subseq(seq_trunc, 1, threshold - 1)
      # shuffle AAStringset to avoid having all the heavy sequences
      # in the last chunk
      seq_set <- sample(c(seq_keep, seq_trunc))
      
      if (all(width(seq_set) < threshold)) return(seq_set)
    }
  }
  
  # ------ Produce encouarging status messages, inputs should be ok.
  # create actual tool name with version number provided
  signalp_version <- paste("signalp", version, sep = '')
  message(paste('Version used...', signalp_version))
  message("running signalp locally...")
  
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
      mp <- suppressMessages(manage_paths(in_path = FALSE,
                             test_mode = signalp_version,
                             path_file = paths))
      full_pa <- mp$path_tibble$path
    } 
    
    # decide how to run signalp depending on a version provided
    
    if (version == 4) {
      # runing signalp versios 4 and 4.1, potentially should work for 5
      
      sp <- tibble::as.tibble(read.table(text = (system(paste(full_pa, "-t",
                                                      organism, out_tmp),
                                                      intern = TRUE))))
      names(sp) <- c("gene_id", "Cmax", "Cpos",
                     "Ymax", "Ypos", "Smax",
                     "Spos", "Smean", "D",
                     "Prediction", "Dmaxcut",
                     "Networks-used")

      # reorder columns to match sp2/3 output:
      sp <- sp %>% select("gene_id", "Cmax", "Cpos",
                          "Ymax", "Ypos", "Smax",
                          "Spos", "Smean", "Prediction")

      sp <- sp %>% filter_( ~ Prediction == 'Y')
      sp$Prediction <- ifelse(sp$Prediction == 'Y', 'Signal peptide')

    } else if (version < 4) {
      # running signalp versions 2 and 3, call parse_signalp for the output
      message('signalp < 4, calling parser for the output...')
      con <- system(paste(full_pa, "-t", organism, out_tmp), intern = TRUE)
      sp <- parse_signalp(input = con, input_type = "system_call")
    }

    message(paste('Candidate sequences with signal peptides...', nrow(sp)))

    if (nrow(sp) == 0) {
      warning('Signal peptide prediction yeilded 0 candidates')
    }

    # generate cropped names for input fasta
    cropped_names <- unname(sapply(names(aaSet), crop_names))
    # replace long names with cropped names
    names(aaSet) <- cropped_names
    # get ids of candidate secreted proteins
    candidate_ids <- sp %>%
                     select_( ~ gene_id) %>%
                     unlist(use.names = FALSE)
    out_fasta_sp <- aaSet[candidate_ids]

    # generate mature sequences
    sp_Cpos <- sp %>% select_( ~ Cpos) %>% unlist(use.names = FALSE)
    cropped_fasta <- subseq(out_fasta_sp, start = sp_Cpos, end = -1)

    # construct output object
    out_obj <- SignalpResult(
      in_fasta = aaSet,
      out_fasta = out_fasta_sp,
      mature_fasta = cropped_fasta,
      sp_version = version,
      sp_tibble = sp)
    
    # check that intended output is valid
    if (validObject(out_obj)) {return(out_obj)}
  } #close simple_fasta
  
  # Handle long sequences if any present in the input,
  # if necessary - run signalp as a parallel process

  # estimate how big is the file, if required - split it into smaller
  # chunks and run signalp as an embarassingly parallel process

  fasta <- truncate_seq(truncate = truncate, fasta, 2000)

  # helper function to estimate approximate length threshold if chunk
  # size exceedes 200000
  estimate_lim <- function(fasta_chunk) {
    len_lim <- (200000 / length(fasta_chunk) + 50)
    if (sum(width(fasta_chunk)) >= 200000) {
      message(paste(
          'fasta size exceedes maximal total residue limit, sequences longer',
          round(len_lim),
          'will be truncated'))
      
      fasta_trunc <- truncate_seq(truncate = truncate, fasta_chunk, len_lim)
      return(fasta_trunc)
    } else {
      return(fasta_chunk)
    }
  }

  if (length(fasta) <= 600) {
    message('Ok for single processing')
    return(simple_signalp(estimate_lim(fasta)))
  } else {
    message('Input fasta contains >600 sequences, entering batch mode...')

    #split and check that chunks do not exceed 200K residue limit
    split_fasta <- sapply(split_XStringSet(fasta, 500), estimate_lim)

    # Initiate cluster
    cl <- makeCluster(cores)
    # run parallel process
    clusterEvalQ(cl, library(SecretSanta))
    clusterExport(cl = cl,
                  varlist = c("paths"),
                  envir = environment())
    result <- parLapply(cl, split_fasta, simple_signalp)
    stopCluster(cl)

    res_comb <- do.call(c, result)
    combined_SignalpResult <- combine_SpResult(unname(res_comb))

    sp_count <- nrow(getSPtibble(combined_SignalpResult))
    message(paste('Candidate sequences with signal peptides...',
                  sp_count))
    if (sp_count == 0) {
      warning('Signal peptide prediction yeilded 0 candidates')
    }

    closeAllConnections()
    return(combined_SignalpResult)
  }
}
