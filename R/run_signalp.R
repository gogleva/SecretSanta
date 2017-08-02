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
#' @export
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
  
  signalp_version <- paste("signalp", version, sep = '')
  message(signalp_version)
  full_pa <- as.character(paths %>% filter(tool == signalp_version) %>% select(path))
    
  # helper function: crop long names for AAStringSet object, return character vector
  crop_names <- function(x){unlist(stringr::str_split(x, " "))[1]}
    
  # ----
  if (version >= 4) {
  # runing signalp versios 4 and 4.1, potentially should work for 5
    sp <- tibble::as.tibble(read.table(text = (system(paste(full_pa, "-t", organism_type, out_tmp), intern = TRUE))))
    names(sp) <- c("gene_id", "Cmax", "Cpos",
                            "Ymax", "Ypos", "Smax",
                            "Spos", "Smean", "D",
                            "Prediction", "Dmaxcut", "Networks-used")
    sp <- sp %>% filter(Prediction == 'Y')
      
  } else if (version < 4) {
  # running signalp versions 2 and 3, call parse_signalp for the output
      message('signalp < 4, calling parser for the output...')
      con <- system(paste(full_pa, "-t", organism_type, out_tmp), intern = TRUE)
      sp <- parse_signalp(input = con, input_type = "system_call")
    }

  # generate cropped names for input fasta
  cropped_names <- unname(sapply(names(fasta), crop_names))
  # replace long names with cropped names
  names(fasta) <- cropped_names
  # get ids of candidate secreted proteins
  candidate_ids <- sp %>% select(gene_id) %>% unlist(use.names = FALSE)
  out_fasta_sp <- fasta[candidate_ids]
    
  # generate mature sequences
    
  sp_Cpos <- sp %>% select(Cpos) %>% unlist(use.names = FALSE)
  cropped_fasta <- subseq(out_fasta_sp, start = sp_Cpos, end = -1)
    
  out_obj <- SignalpResult(in_fasta = fasta,
                           out_fasta = out_fasta_sp, 
                           mature_fasta = cropped_fasta, 
                           sp_version = version,
                           sp_tibble = sp)
  if (validObject(out_obj)) {return(out_obj)}
}
