# HELPER FUNCTIONS FOR PARALLEL SIGNALP:

#' split_XStringSet function
#'
#' This function splits large XStringSet files into chunks of given size and returns a list of AAStringSets, those could be written to tmp files.
#' @param string_set - input AAStringSet that requires spliting;
#' @param chunk_size - number of sequenses in a single chunk;
#' @export
#' @examples 
#' large_aa <- readAAStringSet(system.file("extdata", "Ppalm_prot_ALI_PLTG.fasta", package = "SecretSanta"))
#' split_XStringSet(large_aa, 1000)
#' res <- split_XStringSet(large_aa, 1000)

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

#' combine_SignalpResult function
#'
#' This function combines multiple instances of SignalpResult class, typically generated with parLapply
#' @param arguments - a list of SignalpResult objects to be combined in one
#' @export
#' @examples 
#' inp2 <- CBSResult(in_fasta = readAAStringSet(system.file("extdata", "tail_prot.fasta", package = "SecretSanta")))
#' inp3 <- CBSResult(in_fasta = readAAStringSet(system.file("extdata", "tail2_prot.fasta", package = "SecretSanta")))
#' inp4 <- CBSResult(in_fasta = readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta")))
#' 
#' sp1 <- signalp(input_obj = inp2, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa)
#' sp2 <- signalp(input_obj = inp3, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa)
#' sp3 <- signalp(input_obj = inp4, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa)
#' obj <- list(sp1, sp2, sp3)
#  combined_sp <- combine_SignalpResult(obj)


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
#' @param arguments - a list of CBSResult objects to be combined
#' @export
#' @examples 
#' inp2 <- CBSResult(in_fasta = readAAStringSet(system.file("extdata", "tail_prot.fasta", package = "SecretSanta")))
#' inp3 <- CBSResult(in_fasta = readAAStringSet(system.file("extdata", "tail2_prot.fasta", package = "SecretSanta")))
#' inp4 <- CBSResult(in_fasta = readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta")))
#' combined_CBS <- combine_CBSResult(inp2, inp3, inp4)
#' combine_CBSResult(sp1, sp2, sp3)

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

signalp_parallel <- function(input_obj, version, organism_type, run_mode, paths) {
  
  # ----- Check that inputs are valid
  
  # check that input object belong to CBSResult class
  if (is(input_obj, "CBSResult")) {} else {stop('input_object does not belong to CBSResult superclass')}
  
  # check that supplied runnig mode is valid
  
  if (run_mode %in% c('piper', 'starter')) {} else {stop("Run mode is invalid. Please use 'starter' to initiate prediction pipelie or 'piper' to continue")}
  
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
  
  # Hack - throw away too long seqeunces, longer than 4000 residues
  fasta <- fasta[width(fasta) < 4000]
  warning(paste('Some long sequenses have been thrown away'))
  
  # check signalp versions and organism type
  allowed_versions = c(2,3,4,4.1)
  allowed_organisms = c('euk', 'gram+', 'gram-')
  organism_type <- tolower(organism_type)
  
  signalp_version <- paste("signalp", version, sep = '')
  message(paste('Version used...', signalp_version))
  
  if ((version %in% allowed_versions) & (organism_type %in% allowed_organisms)) {
    message("running signalP locally...")
  } else {
    message('Allowed versions include...:')
    message(cat(allowed_versions))
    message('Allowed organisms iclude...:')
    message(cat(allowed_organisms))
    stop('Input signalp version or specified organism type are invalid.')  
  }
  
  # simplesignalp, takes single AAStringSet as an input and runs signalp on it - function body from run_signalp
  
  simple_signalp <- function(aaSet) { 

    # ---- Run prediction
    # convert fasta to a temporary file:
    out_tmp <- tempfile() #create a temporary file for fasta
    Biostrings::writeXStringSet(aaSet, out_tmp) #write tmp fasta file
    
    # make a system call of signalp based on the tmp file
    
    full_pa <- as.character(paths %>% dplyr::filter(tool == signalp_version) %>% dplyr::select(path))
    
    # helper function: crop long names for AAStringSet object, return character vector
    crop_names <- function(x){unlist(stringr::str_split(x, " "))[1]}
    
    message(paste('Number of submitted sequences...', length(aaSet)))
    
    # ----
    if (version >= 4) {
      # runing signalp versios 4 and 4.1, potentially should work for 5
      sp <- tibble::as.tibble(read.table(text = (system(paste(full_pa, "-t", organism_type, out_tmp), intern = TRUE))))
      names(sp) <- c("gene_id", "Cmax", "Cpos",
                     "Ymax", "Ypos", "Smax",
                     "Spos", "Smean", "D",
                     "Prediction", "Dmaxcut", "Networks-used")
      # reorder columns to match sp2/3 output:
      
      sp <- sp %>% dplyr::select("gene_id",
                                 "Cmax",
                                 "Cpos",
                                 "Ymax",
                                 "Ypos",
                                 "Smax",
                                 "Spos",
                                 "Smean",
                                 "Prediction")
      
      sp <- sp %>% dplyr::filter(Prediction == 'Y')
      sp <- dplyr::mutate(sp, Prediction = ifelse(Prediction == 'Y', 'Signal peptide'))
      
      
    } else if (version < 4) {
      # running signalp versions 2 and 3, call parse_signalp for the output
      message('signalp < 4, calling parser for the output...')
      con <- system(paste(full_pa, "-t", organism_type, out_tmp), intern = TRUE)
      sp <- parse_signalp(input = con, input_type = "system_call")
    }
    
    message(paste('Number of candidate sequences with signal peptides...', nrow(sp)))
    
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
    
    # costruct output object
    
    out_obj <- SignalpResult(in_fasta = aaSet,
                             out_fasta = out_fasta_sp,
                             mature_fasta = cropped_fasta,
                             sp_version = version,
                             sp_tibble = sp)
    if (validObject(out_obj)) {return(out_obj)}
  }
  
  # estimate how big is the file, if required - split it into smaller chunks and run
  # signalp as an embarassingly parallel process
  
  if (length(fasta) < 500) {message('Ok for single processing')
    simple_signalp(fasta)
  } else {
    message('Input fasta contains >500 sequences, entering batch mode...')
    split_fasta <- split_XStringSet(fasta, 500)
    
    # Calculate the number of cores
    no_cores <- detectCores()
   
    # Initiate cluster
    cl <- makeCluster(no_cores)
    # run parallel process
    
    clusterEvalQ(cl, library("SecretSanta"))
    clusterExport(cl=cl, varlist=c("my_pa")) # or path?
    result <- parLapply(cl, split_fasta, simple_signalp)
  
    stopCluster(cl)
    res_comb <- do.call(c,result)
    return(combine_SignalpResult(unname(res_comb)))
  }
}
  

# to do: need to clean tmp files on exit signalp_chunk

# test run:
 
my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))

aa_1K <- readAAStringSet("/home/anna/anna/Labjournal/SecretSanta_external/test_fastas/medium_1K.fasta")
inp_1K <- CBSResult(in_fasta = aa_1K)

aa_2K <- readAAStringSet("/home/anna/anna/Labjournal/SecretSanta_external/test_fastas/medium_2K.fasta")
inp_2K <- CBSResult(in_fasta = aa_2K)

sp_1K_par <- signalp_parallel(inp_1K, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa)

# profile 1K input, parallel
microbenchmark::microbenchmark(signalp_parallel(inp_1K, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa), times = 1)
microbenchmark::microbenchmark(signalp_parallel(inp_1K, version = 3, organism_type = 'euk', run_mode = 'starter', paths = my_pa), times = 1)
microbenchmark::microbenchmark(signalp_parallel(inp_1K, version = 4, organism_type = 'euk', run_mode = 'starter', paths = my_pa), times = 1)

# profile 2K input, parallel
microbenchmark::microbenchmark(signalp_parallel(inp_2K, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa), times = 1)
microbenchmark::microbenchmark(signalp_parallel(inp_2K, version = 3, organism_type = 'euk', run_mode = 'starter', paths = my_pa), times = 1)
microbenchmark::microbenchmark(signalp_parallel(inp_2K, version = 4, organism_type = 'euk', run_mode = 'starter', paths = my_pa), times = 1)


# profile 2K input, non-parallel
microbenchmark::microbenchmark(signalp(inp_2K, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa), times = 1) # fails due to the input limits
microbenchmark::microbenchmark(signalp(inp_2K, version = 3, organism_type = 'euk', run_mode = 'starter', paths = my_pa), times = 1)
microbenchmark::microbenchmark(signalp(inp_2K, version = 4, organism_type = 'euk', run_mode = 'starter', paths = my_pa), times = 1)

# profile 1K input, non-parallel
microbenchmark::microbenchmark(signalp(inp_2K, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa), times = 1) # fails due to the input limits
microbenchmark::microbenchmark(signalp(inp_2K, version = 3, organism_type = 'euk', run_mode = 'starter', paths = my_pa), times = 1)
microbenchmark::microbenchmark(signalp(inp_2K, version = 4, organism_type = 'euk', run_mode = 'starter', paths = my_pa), times = 1)
