# HELPER FUNCTIONS FOR PARALLEL SIGNALP:

#' split_XStringSet function
#'
#' This function splits large XStringSet files into chunks of given size and returns a list of AAStringSets, those could be written to tmp files.
#' @param string_set - input AAStringSet that requires spliting;
#' @param chunk_size - number of sequenses in a single chunk;
#' @export

split_XStringSet <- function(string_set, chunk_size){
                                
                                if (!(is(string_set, 'XStringSet'))) {
                                  stop('Input string_set does not belong to XStringSet class')
                                  }
                                
                                lst <- length(string_set)
                                
                                if (chunk_size > lst) {
                                  stop('Chunk size exceeds total seq number')
                                }
                                
                                total_seq  <- c(1:lst)
                                chunks <- split(total_seq, ceiling(seq_along(total_seq)/chunk_size))
                                seq_chunker <- function(x) {chunk <- string_set[x]}
                                lapply(chunks, seq_chunker) 
                                
}

#combine_SignalpResult function
#'
#' This function combines multiple instances of SignalpResult class, typically generated with parLapply
#' @param arguments - a list of SignalpResult objects to be combined in one
#' @export
#' @examples 
#' my_pa <- manage_paths(system.file(
#'                       "extdata",
#'                       "sample_paths",
#'                        package = "SecretSanta"))
#'                                                  
#' inp2 <- CBSResult(in_fasta =
#'                   readAAStringSet(
#'                   system.file("extdata",
#'                   "tail_prot.fasta",
#'                   package = "SecretSanta")
#'                   ))
#' inp3 <- CBSResult(in_fasta = 
#'                   readAAStringSet(
#'                   system.file("extdata",
#'                   "tail2_prot.fasta",
#'                   package = "SecretSanta")
#'                   ))
#' inp4 <- CBSResult(in_fasta =
#'                  readAAStringSet(
#'                  system.file("extdata",
#'                  "sample_prot_100.fasta",
#'                  package = "SecretSanta")
#'                  ))
#'                   
#' sp1 <- signalp(input_obj = inp2,
#'                version = 2,
#'                organism_type = 'euk',
#'                run_mode = 'starter',
#'                paths = my_pa)
#' sp2 <- signalp(input_obj = inp3,
#'                version = 2,
#'                organism_type = 'euk',
#'                run_mode = 'starter',
#'                paths = my_pa)
#' sp3 <- signalp(input_obj = inp4,
#'                version = 2,
#'                organism_type = 'euk',
#'                run_mode = 'starter', 
#'                paths = my_pa)
#' obj <- list(sp1, sp2, sp3)
#' combined_sp <- combine_SignalpResult(obj)


combine_SignalpResult <- function(arguments) {
                                    if (all(sapply(arguments, is, 'SignalpResult'))) {
                                    } else {                               
                                      stop('Some objects from arguments list do not belong to SignalpResult class.')
                                    }
                                    
                                    c_in_fasta <- do.call(c, (lapply(arguments, getInfasta)))
                                    c_out_fasta <- do.call(c, (lapply(arguments, getOutfasta)))
                                    c_mature_fasta <- do.call(c, (lapply(arguments, getMatfasta)))
                                    c_sp_tibble <- do.call(rbind, (lapply(arguments, getSPtibble)))
                                    c_sp_version <- unlist((lapply(arguments, getSPversion))[1])

                                    c_obj <- SignalpResult(in_fasta = c_in_fasta,
                                                       out_fasta = c_out_fasta,
                                                       mature_fasta = c_mature_fasta,
                                                       sp_tibble = c_sp_tibble,
                                                       sp_version = c_sp_version)
}

#' combine_CBSResult function
#'
#' Minimal function to combine objects of CBSResult class.
#' @param ... - a list of CBSResult objects to be combined
#' @export
#' @examples 
#' inp2 <- CBSResult(in_fasta = 
#'                   readAAStringSet(
#'                   system.file("extdata",
#'                   "tail_prot.fasta",
#'                   package = "SecretSanta")
#'                   ))
#' inp3 <- CBSResult(in_fasta =
#'                   readAAStringSet(
#'                   system.file("extdata",
#'                   "tail2_prot.fasta",
#'                   package = "SecretSanta")
#'                   ))
#' inp4 <- CBSResult(in_fasta =
#'                   readAAStringSet(
#'                   system.file("extdata",
#'                   "sample_prot_100.fasta",
#'                   package = "SecretSanta")))
#' combined_CBS <- combine_CBSResult(inp2, inp3, inp4)

combine_CBSResult <- function(...) {
                              arguments <- list(...)
                              
                              if (all(sapply(arguments, is, 'CBSResult'))) {
                              } else {                               
                                stop('Some or all objects from the arguments list do not belong to CBSResult class.')
                              }
                              
                              if (any(sapply(arguments, is, 'SignalpResult'))) {
                                warning('Only in_fasta and out_fasta slots will be combined')
                              }
                              
                              comb_in_fasta <- do.call(c, (lapply(arguments, getInfasta)))
                              comb_out_fasta <- do.call(c, (lapply(arguments, getOutfasta)))
                              c_obj <- CBSResult(in_fasta = comb_in_fasta,
                                                 out_fasta = comb_out_fasta)

}

# parallel version of signalp:

#' signalp function
#'
#' This function calls local signalp to predict the presence and location of signal peptide cleavage sites in amino acid sequences; automatically splits large input files (>500 sequnces) and runs signalp prediction as an embarassingly parallel process on all the CPUs available.
#' @param input_obj   an instance of CBSResult class containing protein sequences as on of the attributes
#' @param version  signalp version to run, supported versions:
#' \itemize{
#' \item 2
#' \item 3
#' \item 4
#' \item 4.1
#' } 
#' @param organism_type possible values: 
#' \itemize{
#' \item euk - for eukaryotes;
#' \item gram+ - for gram-positive bacteria;
#' \item gram- - for gram-negative bacteria;
#' }
#' @param run_mode
#' \itemize{
#' \item starter - if it is the first step in pipeline;
#' \item piper - if you run this function on the output of other CBS tools;
#' }
#' @param paths   tibble with paths to external dependencies, generated with \code{\link{manage_paths}} function
#' @param truncate logical, if TRUE - sequences longer 2000 residues will be truncated to this length limit and renamed. If FALSE - long sequences will be excluded from the analysis. Default = TRUE.
#' @return an object of SignalpResult class
#' @export
#' @examples
#' 
#' #Example pipe would loook like this:
#' 
#' # set paths for external dependencies with manage_paths()
#' my_pa <- manage_paths(system.file("extdata",
#'                                   "sample_paths",
#'                                    package = "SecretSanta"))
#' 
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata",
#'                                   "sample_prot_100.fasta",
#'                                    package = "SecretSanta"))
#' 
#' # assign this object to the input_fasta slot of empty CBSResult object
#' inp <- CBSResult(in_fasta = aa[1:10])
#' 
#' # run signalp2 on the initial file:
#' step1_sp2 <- signalp(inp,
#'                      version = 2,
#'                      'euk', 
#'                      run_mode = "starter",
#'                      paths = my_pa)
#' 
#' # run signalp3 on the result object, will automatically pass out_fasta slot to signalp3:
#' step2_sp3 <- signalp(step1_sp2,
#'                      version = 3,
#'                      'euk',
#'                      run_mode = "piper",
#'                      paths = my_pa)
#' 
#' # run signalp4 on the result object, will automatically pass out_fasta slot to signalp4:
#' step3_sp4 <- signalp(step2_sp3,
#'                      version = 4,
#'                      'euk',
#'                      run_mode = "piper",
#'                      paths = my_pa)

signalp <- function(input_obj, version, organism_type, run_mode, paths, truncate = NULL) {
  
  # ----- Set default value for parameters if not provided:
  
  if (is.null(truncate)) truncate = TRUE else truncate
  if (is.logical(truncate)) {} else {stop('truncate parameter must be logical')}
  
  # ------ Helper functions:
  
  # helper function: crop long names for AAStringSet object, return character vector
  crop_names <- function(x){unlist(stringr::str_split(x, " "))[1]}
  
  # helper function to truncate log sequences or throw them away

  truncate_seq <- function(truncate, seq_set, threshold) {
    drop_n <- length(seq_set[width(seq_set) >= threshold])
    
    if (drop_n == 0) return(seq_set) 

    if (truncate == F) {
      seq_set <- seq_set[width(seq_set) < threshold]
      warning(paste(drop_n, 'long sequenses have been thrown away'))
      return(seq_set)

    } else if (truncate == T) {
      message(paste(drop_n, 'sequences to be truncated'))
      seq_keep <- seq_set[width(seq_set) < threshold] # not so long sequences
      seq_trunc <- seq_set[width(seq_set) >= threshold] # sequences we need to truncate
      t_names <- paste(unname(sapply(names(seq_trunc), crop_names)),
                       '_truncated',
                       sep = '')
      names(seq_trunc) <- t_names #new names for sequences to be truncated
      seq_trunc <- Biostrings::subseq(seq_trunc, 1, threshold - 1)
      seq_set <- sample(c(seq_keep, seq_trunc)) # shuffle AAStringset to avoid having all the heavy sequences in the last chunk

      if (all(width(seq_set) < threshold)) return(seq_set)
    }
  }
  
  # ----- Check that inputs are valid

  # check that input object belong to CBSResult class
  if (is(input_obj, "CBSResult")) {} else {stop('input_object does not belong to CBSResult superclass')}

  # check that supplied runnig mode is valid

  if (run_mode %in% c('piper', 'starter')) {} else {
    stop("Run mode is invalid. Please use 'starter' to initiate prediction pipelie or 'piper' to continue")}

  # check that input_object contains non-empty in/out_fasta for starter/piper

  if (run_mode == 'starter') {
    if (length(getInfasta(input_obj)) != 0) {
      fasta <- getInfasta(input_obj)
    } else {stop('in_fasta attribute is empty')}
  } else if (run_mode == 'piper') {
    if (length(getOutfasta(input_obj)) != 0) {
      fasta <- getOutfasta(input_obj)
    } else {stop('out_fasta attribute is empty')}
  }

  # check signalp versions and organism type
  allowed_versions = c(2,3,4,4.1)
  allowed_organisms = c('euk', 'gram+', 'gram-')
  organism_type <- tolower(organism_type)
  
  # 
  signalp_version <- paste("signalp", version, sep = '')
  message(paste('Version used...', signalp_version))

  if ((version %in% allowed_versions) & (organism_type %in% allowed_organisms)) {
    message("running signalp locally...")
  } else {
    stop('Input signalp version or specified organism type are invalid.')
  }

  # simple signalp, takes single AAStringSet as an input and runs signalp on it - function body from run_signalp

  simple_signalp <- function(aaSet) {

    # ---- Run prediction
    # convert fasta to a temporary file:
    out_tmp <- tempfile() #create a temporary file for fasta
    Biostrings::writeXStringSet(aaSet, out_tmp) #write tmp fasta file

    # make a system call of signalp based on the tmp file

    full_pa <- as.character(paths %>% dplyr::filter(tool == signalp_version) %>% dplyr::select(path))
    message(paste('Submitted sequences...', length(aaSet)))

    # ----
    if (version >= 4) {
      # runing signalp versios 4 and 4.1, potentially should work for 5
      sp <- tibble::as.tibble(read.table(text = (system(paste(full_pa, "-t", organism_type, out_tmp), intern = TRUE))))
      names(sp) <- c("gene_id", "Cmax", "Cpos",
                     "Ymax", "Ypos", "Smax",
                     "Spos", "Smean", "D",
                     "Prediction", "Dmaxcut", "Networks-used")
      # reorder columns to match sp2/3 output:

      sp <- sp %>% dplyr::select("gene_id", "Cmax", "Cpos",
                                 "Ymax", "Ypos", "Smax",
                                 "Spos", "Smean", "Prediction")

      sp <- sp %>% dplyr::filter(Prediction == 'Y')
      sp <- dplyr::mutate(sp, Prediction = ifelse(Prediction == 'Y', 'Signal peptide'))


    } else if (version < 4) {
      # running signalp versions 2 and 3, call parse_signalp for the output
      message('signalp < 4, calling parser for the output...')
      con <- system(paste(full_pa, "-t", organism_type, out_tmp), intern = TRUE)
      sp <- parse_signalp(input = con, input_type = "system_call")
    }

    message(paste('Candidate sequences with signal peptides...', nrow(sp)))

    if (nrow(sp) == 0) {warning('Signal peptide prediction yeilded 0 candidates')}

    # generate cropped names for input fasta
    cropped_names <- unname(sapply(names(aaSet), crop_names))
    # replace long names with cropped names
    names(aaSet) <- cropped_names
    # get ids of candidate secreted proteins
    candidate_ids <- sp %>% dplyr::select(gene_id) %>% unlist(use.names = FALSE)
    out_fasta_sp <- aaSet[candidate_ids]

    # generate mature sequences

    sp_Cpos <- sp %>% dplyr::select(Cpos) %>% unlist(use.names = FALSE)
    cropped_fasta <- subseq(out_fasta_sp, start = sp_Cpos, end = -1)

    # construct output object

    out_obj <- SignalpResult(in_fasta = aaSet,
                             out_fasta = out_fasta_sp,
                             mature_fasta = cropped_fasta,
                             sp_version = version,
                             sp_tibble = sp)
    if (validObject(out_obj)) {return(out_obj)}
  }

  # estimate how big is the file, if required - split it into smaller chunks and run
  # signalp as an embarassingly parallel process

  # Handle long sequences if any present in the input, id necessary - run signalp as 
  # an embarassingly parallel process
  
  fasta <- truncate_seq(truncate = truncate, fasta, 2000)

  # to do: check total number of residues
  
  # helper function to estimate approximate length threshold if chink size exceedes 200000
  estimate_lim <- function(fasta_chunk){
    len_lim <- (200000/ length(fasta_chunk) + 50)
    if (sum(width(fasta_chunk)) >= 200000) {
      message(paste('fasta size exceedes maximal total residue limit, sequences longer',
                    round(len_lim),
                    'will be truncated'))
      fasta_trunc <- truncate_seq(truncate = truncate, fasta_chunk, len_lim)
      return(fasta_trunc)
    } else {
      return(fasta_chunk)
    }
  }
  
  if (length(fasta) <= 600) {message('Ok for single processing')
    return(simple_signalp(estimate_lim(fasta)))
  } else {
    message('Input fasta contains >600 sequences, entering batch mode...')
    split_fasta <- sapply(split_XStringSet(fasta, 500), estimate_lim) #split and check that chunks do not exceed 200K residue limit

    # Calculate the number of cores
    no_cores <- detectCores()

    # Initiate cluster
    cl <- makeCluster(no_cores)
    # run parallel process

    clusterEvalQ(cl, library("SecretSanta"))
    clusterExport(cl=cl, varlist=c("paths"), envir = environment()) #
    result <- parLapply(cl, split_fasta, simple_signalp)
    stopCluster(cl)

    res_comb <- do.call(c,result)
    combined_SignalpResult <- combine_SignalpResult(unname(res_comb))
    
    sp_count <- nrow(getSPtibble(combined_SignalpResult))    
    message(paste('Candidate sequences with signal peptides...', sp_count))
    if (sp_count == 0) {warning('Signal peptide prediction yeilded 0 candidates')}
    
    closeAllConnections()
    return(combined_SignalpResult)
  }
   
}
