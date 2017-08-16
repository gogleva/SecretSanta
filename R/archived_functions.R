#' Parse_signalP_web function
#'
#' This function scraps web output of SignalP tool and returns output in short format organised in a dataframe
#' @param url url used to produce the output
#' Parse_signalP_web()

Parse_signalP_web <- function(url){
  webpage <- xml2::read_html(url)
  raw_data_html <- rvest::html_nodes(webpage, 'pre')
  data <- rvest::html_text(raw_data_html)
  raw <- unlist(strsplit(data, '\n'))
  le <- length(raw)
  d <- raw[9:le - 4] #draft parameters
  d2 <- as.data.frame(do.call("rbind", d2)) 
  d2
}

#' pipe_signalp function
#'
#' This function pipes outputs from different version of signalp
#' @param piping_sequence  specify piping sequence: vector of signalp versions \cr
#'                         accepted formats: c(1,2,3) or c("signalp2", "signalp3", "signalp4")
#' @param input_fasta fasta file with proteins to initiate the pipe
#' @examples
#' piping_seq <- c("signalp2", "signalp3", "signalp4")
#' piping_seq_num <- c(2,3,3)
#' pipe_signalp(piping_sequence = piping_seq, 'some_fasta')
#' pipe_signalp(piping_sequence = piping_seq_num, 'some_fasta')
#'

pipe_signalp <- function(piping_sequence, input_fasta) {
  # check that pipong sequence does not contain duplicates
  if (length(unique(piping_sequence)) == length(piping_sequence)) {
    # generate starting message
    message("Hi, let's pipe!")
    if (all(stringr::str_detect(piping_sequence, "signalp"))) {
      piping_sequence
    } else if(is.numeric(piping_sequence)){
      piping_sequence <- paste('signalp', piping_seq_num, sep = '')
    }
    message(paste(piping_sequence[-length(piping_sequence)], '--> '), piping_sequence[length(piping_sequence)])
    message('pipeline includes ', length(piping_seq), ' steps')
  }else{
    message('Aborting..')
    stop('please make sure that there are no duplicte tools in the specified piping sequence')
  }  
  
  # piper:
  
  # run signalp#1(input_fasta) -> outputs tibble -> generate (fasta)'
  # run signalp#2(fasta') -> output tibble' -> generate (fasta)''
  # run signalp#3(fasta'') -> output tibble'' == result
  
  
}

#helper recursive piper function: 

#piper <- function(some_sequence) {
#   if  
#}


#to do: teach run_signalp to output fasta file and/or tibble with results

#' Convert \code{data.frame} to \code{list}.
#' 
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @param x A \code{data.frame} object.
#' @examples
#' my_result <- foo(iris)
#'
foo <- function(x) {
  x %>%
    as.list()
}



#' signalp function
#'
#' This function calls local signalp to predict the presence and location of signal peptide cleavage sites in amino acid sequences.
#' @param input_object    an instance of CBSResult class containing protein sequences as on of the attributes
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
#' @return an object of SignalpResult class
#' 
#' @examples
#' Example pipe would loook like this:
#' 
#' # set paths for external dependencies with manage_paths()
#' my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
#' 
#' # initialise SignalpResult object
#' inp <- SignalpResult()
#' 
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), use.names = TRUE)
#' 
#' # assign this object to the input_fasta slot of SignalpResult object
#' inp <- setInfasta(inp, aa)
#' 
#' # run signalp2 on the initial file:
#' step1_sp2 <- signalp(inp, version = 2, 'euk', run_mode = "starter", paths = my_pa)
#' 
#' # run signalp3 on the result object, will automatically pass out_fasta slot to signalp3:
#' step2_sp3 <- signalp(step1_sp2, version = 3, 'euk', run_mode = "piper", paths = my_pa)
#' 
#' # run signalp4 on the result object, will automatically pass out_fasta slot to signalp4:
#' step3_sp4 <- signalp(step2_sp3, version = 4, 'euk', run_mode = "piper", paths = my_pa)


signalp <- function(input_obj, version, organism_type, run_mode, paths) {
  
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
  
  # ---- Run prediction
  # convert fasta to a temporary file:
  out_tmp <- tempfile() #create a temporary file for fasta
  Biostrings::writeXStringSet(fasta, out_tmp) #write tmp fasta file
  
  # make a system call of signalp based on the tmp file
  
  full_pa <- as.character(paths %>% dplyr::filter(tool == signalp_version) %>% dplyr::select(path))
  
  # helper function: crop long names for AAStringSet object, return character vector
  crop_names <- function(x){unlist(stringr::str_split(x, " "))[1]}
  
  message(paste('Number of submitted sequences...', length(fasta)))
  
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
  cropped_names <- unname(sapply(names(fasta), crop_names))
  # replace long names with cropped names
  names(fasta) <- cropped_names
  # get ids of candidate secreted proteins
  candidate_ids <- sp %>% dplyr::select(gene_id) %>% unlist(use.names = FALSE)
  out_fasta_sp <- fasta[candidate_ids]
  
  # generate mature sequences
  
  sp_Cpos <- sp %>% dplyr::select(Cpos) %>% unlist(use.names = FALSE)
  cropped_fasta <- subseq(out_fasta_sp, start = sp_Cpos, end = -1)
  
  # costruct output object
  
  out_obj <- SignalpResult(in_fasta = fasta,
                           out_fasta = out_fasta_sp, 
                           mature_fasta = cropped_fasta, 
                           sp_version = version,
                           sp_tibble = sp)
  if (validObject(out_obj)) {return(out_obj)}
}



#' targetp function
#'
#' This function calls local targetp to predict subcellular localisation of a protein.
#' @param input_object    an instance of CBSResult class containing protein sequences as on of the attributes
#' @param network_type possible values: 
#' \itemize{
#' \item P - for plants;
#' \item N - for non-plants;
#' }
#' @param run_mode
#' \itemize{
#' \item starter - if it is the first step in pipeline;
#' \item piper - if you run this function on the output of other CBS tools;
#' }
#' @param paths   tibble with paths to external dependencies, generated with \code{\link{manage_paths}} function
#' @return an object of TargetpResult class
#' @examples 
#' my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
#' # initialise SignalpResult object
#' inp <- SignalpResult()
#' 
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), use.names = TRUE)
# 
#' # assign this object to the input_fasta slot of SignalpResult object
#' inp <- setInfasta(inp, aa)
#' 
#' tp_result <- targetp(input_object = inp, network_type = 'N', run_mode = 'starter', paths = my_pa)

targetp <- function(input_object, network_type, run_mode, paths) {
  # ----- Check that inputs are valid
  
  # check that input object belong to CBSResult class
  
  if (is(input_object, "CBSResult")) {} else {stop('input_object does not belong to CBSResult superclass')}
  
  # check that supplied runnig mode is valid
  
  if (run_mode %in% c('piper', 'starter')) {} else {stop("Run mode is invalid. Please use 'starter' to initiate prediction pipelie or 'piper' to continue")}
  
  # check that input_object contains non-empty in/out_fasta for starter/piper
  
  if (run_mode == 'starter') {
    if (length(getInfasta(input_object)) != 0) {
      fasta <- getInfasta(input_object)
    } else {stop('in_fasta attribute is empty')}
  } else if (run_mode == 'piper') {
    if (length(getOutfasta(input_object)) != 0) {
      fasta <- getOutfasta(input_object)
    } else {stop('out_fasta attribute is empty')}
  }
  
  allowed_networks = c('P', 'N')
  
  if (network_type %in% allowed_networks) {
    message("running targetp locally...")
  } else {
    stop('Specified network_type is invalid.')  
  }
  
  #----- Run targetp prediction:
  
  # convert fasta to a temporary file:
  out_tmp <- tempfile() #create a temporary file for fasta
  Biostrings::writeXStringSet(fasta, out_tmp) #write tmp fasta file
  
  # make a system call of targetp based on the tmp file
  
  full_pa <- as.character(paths %>% dplyr::filter(tool == 'targetp') %>% dplyr::select(path))
  
  # helper function: crop long names for AAStringSet object, return character vector
  crop_names <- function(x){unlist(stringr::str_split(x, " "))[1]}
  message(paste('Number of submitted sequences...', length(fasta)))
  
  # generate cropped names for input fasta
  cropped_names <- unname(sapply(names(fasta), crop_names))
  # replace long names with cropped names
  names(fasta) <- cropped_names
  
  #run targetp:
  
  NN <- paste('-', network_type, sep = '')
  tp <- tibble::as.tibble(read.table(text = (system(paste(full_pa, NN, out_tmp), intern = TRUE)[1: length(fasta) + 8])))
  if (network_type == 'N') {
    names(tp) <- c('gene_id', 'length', 'mTP', 'sp', 'other', 'TP_localization', 'RC')}
  else if (network_type == 'P') {
    names(tp) <- c('gene_id', 'length', 'cTP', 'mTP', 'sp', 'other', 'TP_localization', 'RC')  
  }
  
  tp <- tp %>% dplyr::filter(TP_localization == 'S')
  
  message(paste('Number of candidate secreted sequences', nrow(tp)))         
  
  candidate_ids <- tp %>% dplyr::select(gene_id) %>% unlist(use.names = FALSE)
  out_fasta_tp <- fasta[candidate_ids]
  
  # generate output object:
  
  out_obj <- TargetpResult(in_fasta = fasta,
                           out_fasta = out_fasta_tp,
                           tp_tibble = tp
  )
  if (validObject(out_obj)) {return(out_obj)}
}

