#' split_XStringSet function
#'
#' This function splits large XStringSet files into chunks of given size and writes them in tmp files
#' @param string_set
#' @param chunk_size
#' @param prefix  file name prefix for tmp files
#' @export
#' @examples 
#' large_aa <- readAAStringSet(system.file("extdata", "Ppalm_prot_ALI_PLTG.fasta", package = "SecretSanta"))
#' split_XStringSet(large_aa, 1000, 'test')

split_XStringSet <- function(string_set, chunk_size, prefix){
  
  total_seq  <- c(1:length(string_set))
  chunks <- split(total_seq, ceiling(seq_along(total_seq)/chunk_size))
  seq_chunker <- function(x) {
                              chunk <- string_set[x]
                              out_tmp <- tempfile(pattern = prefix)
                              Biostrings::writeXStringSet(chunk, out_tmp)
  }
  
  invisible(lapply(chunks, seq_chunker))
}




# ### experiments with parallelisation
# 
# library(parallel)
# 
# # function to split tmp fasta into multiple fastas
# # run signalp for each of them
# # combine the result
# # lapply is my friend
# 
# # Calculate the number of cores
# no_cores <- detectCores() - 1
# # Initiate cluster
# cl <- makeCluster(no_cores)
# # large real-life file:
# 
# parLapply()
# stopCluster(cl)


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
  
  # estimate how big is the file, if required - split it into smaller chunks
  
  if (length(fasta) < 500) {message('Ok for single processing')
    } else {
      message('Input fasta is to big, entering batch mode...')
    }
  
}

large_aa <- readAAStringSet(system.file("extdata", "Ppalm_prot_ALI_PLTG.fasta", package = "SecretSanta"))
large_inp <- CBSResult(in_fasta = readAAStringSet(system.file("extdata", "Ppalm_prot_ALI_PLTG.fasta", package = "SecretSanta")))

signalp_parallel(large_inp, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa)

  # 
  # 
  # 
  # # ---- Run prediction
  # # convert fasta to a temporary file:
  # out_tmp <- tempfile() #create a temporary file for fasta
  # Biostrings::writeXStringSet(fasta, out_tmp) #write tmp fasta file
  # 
  # # make a system call of signalp based on the tmp file
  # 
  # full_pa <- as.character(paths %>% dplyr::filter(tool == signalp_version) %>% dplyr::select(path))
  # 
  # # helper function: crop long names for AAStringSet object, return character vector
  # crop_names <- function(x){unlist(stringr::str_split(x, " "))[1]}
  # 
  # message(paste('Number of submitted sequences...', length(fasta)))
  # 
  # # ----
  # if (version >= 4) {
  #   # runing signalp versios 4 and 4.1, potentially should work for 5
  #   sp <- tibble::as.tibble(read.table(text = (system(paste(full_pa, "-t", organism_type, out_tmp), intern = TRUE))))
  #   names(sp) <- c("gene_id", "Cmax", "Cpos",
  #                  "Ymax", "Ypos", "Smax",
  #                  "Spos", "Smean", "D",
  #                  "Prediction", "Dmaxcut", "Networks-used")
  #   # reorder columns to match sp2/3 output:
  #   
  #   sp <- sp %>% dplyr::select("gene_id",
  #                              "Cmax",
  #                              "Cpos",
  #                              "Ymax",
  #                              "Ypos",
  #                              "Smax",
  #                              "Spos",
  #                              "Smean",
  #                              "Prediction")
  #   
  #   sp <- sp %>% dplyr::filter(Prediction == 'Y')
  #   sp <- dplyr::mutate(sp, Prediction = ifelse(Prediction == 'Y', 'Signal peptide'))
  #   
  #   
  # } else if (version < 4) {
  #   # running signalp versions 2 and 3, call parse_signalp for the output
  #   message('signalp < 4, calling parser for the output...')
  #   con <- system(paste(full_pa, "-t", organism_type, out_tmp), intern = TRUE)
  #   sp <- parse_signalp(input = con, input_type = "system_call")
  # }
  # 
  # message(paste('Number of candidate sequences with signal peptides...', nrow(sp)))
  # 
  # if (nrow(sp) == 0) {warning('Signal peptide prediction yeilded 0 candidates')}
  # 
  # # generate cropped names for input fasta
  # cropped_names <- unname(sapply(names(fasta), crop_names))
  # # replace long names with cropped names
  # names(fasta) <- cropped_names
  # # get ids of candidate secreted proteins
  # candidate_ids <- sp %>% dplyr::select(gene_id) %>% unlist(use.names = FALSE)
  # out_fasta_sp <- fasta[candidate_ids]
  # 
  # # generate mature sequences
  # 
  # sp_Cpos <- sp %>% dplyr::select(Cpos) %>% unlist(use.names = FALSE)
  # cropped_fasta <- subseq(out_fasta_sp, start = sp_Cpos, end = -1)
  # 
  # # costruct output object
  # 
  # out_obj <- SignalpResult(in_fasta = fasta,
  #                          out_fasta = out_fasta_sp, 
  #                          mature_fasta = cropped_fasta, 
  #                          sp_version = version,
  #                          sp_tibble = sp)
  # if (validObject(out_obj)) {return(out_obj)}
#}


# need to write a function that splits XstringSet into cmaller chunks of given length


large_aa <- readAAStringSet(system.file("extdata", "Ppalm_prot_ALI_PLTG.fasta", package = "SecretSanta"))


