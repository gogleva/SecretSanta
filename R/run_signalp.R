#' signalp function
#'
#' This function calls local SignalP
#' Please ensure that respective version of SignalP is downloaded, installed and respective path is added to $PATH variable
#' @param input    input object (any from SecretSanta objects family) with protein sequences as on of the attributes
#' @param version  signalp version to run, allowed  values: 2, 3, 4, 4.1 
#' @param organism_type Allowed values: 'euk', 'gram+', 'gram-'                 
#' @export
#' @examples
#' result <- signalp(proteins = "SecretSanta/inst/extdata/sample_prot.fasta", organism_type = 'euk', version = 4) #old
#' w <- SignalpResult()
#' class(w)
#' w <- setInfasta(w, bb) #instance of the SignalpResult object with populated Infasta attribute
#' getInfasta(w)
#' r4 <- signalp(w, version = 4, 'euk')
#' r3 <- signalp(w, version = 3, 'euk')
#' r2 <- signalp(w, version = 2, 'euk')
#' r2 <- signalp(w, version = 4, organism_type = 'gra')

signalp <- function(input_obj, version, organism_type) {
  # validity checks for signalp version and organism_type inputs
  allowed_versions = c(2,3,4,4.1)
  allowed_organisms = c('euk', 'gram+', 'gram-')
  
  if ((version %in% allowed_versions) & (organism_type %in% allowed_organisms)){
    message("running signalP locally...")
    # get fasta from SignalpResult object
    fasta <- getInfasta(input_obj)
    # convert it to a temporary file:
    out_tmp <- tempfile() #create a temporary file for fasta
    Biostrings::writeXStringSet(fasta, out_tmp) #write tmp fasta file
    # make a system call of signalp based on the tmp file
    signalp_version <- paste("signalp", version, sep = '')
    message(signalp_version)
    full_pa <- as.character(secret_paths %>% filter(tool == signalp_version) %>% select(path))
    
    # helper function: crop long names for AAStringSet object, return character vector
    crop_names <- function(x){unlist(stringr::str_split(x, " "))[1]}
    
    # ----
    if (version >= 4) {
    # runing signalp versios 4 and 4.1    
      sp <- tibble::as.tibble(read.table(text = (system(paste(full_pa, "-t", organism_type, out_tmp), intern = TRUE))))
      names(sp) <- c("gene_id", "Cmax", "Cpos",
                            "Ymax", "Ypos", "Smax",
                            "Spos", "Smean", "D",
                            "Prediction", "Dmaxcut", "Networks-used")
      sp <- sp %>% filter(Prediction == 'Y')
      
    } else if (version < 4) {
    # running signalp versions 2 and 3, call parser for the initial output
      message('signalp < 4, calling parser for the output...')
      con <- system(paste(full_pa, "-t", organism_type, out_tmp), intern = TRUE)
      sp <- parse_signalp(input = con, input_type = "system_call")
    }

    #generate cropped names for input fasta
    cropped_names <- unname(sapply(names(fasta), crop_names))
    #replace long names with cropped names
    names(fasta) <- cropped_names
    #get ids of candidate secreted proteins
    candidate_ids <- sp %>% select(gene_id) %>% unlist(use.names = FALSE)
    out_fasta_sp <- fasta[candidate_ids]
    
    #generate mature sequences
    
    sp_Cpos <- sp %>% select(Cpos) %>% unlist(use.names = FALSE)
    cropped_fasta <- subseq(out_fasta_sp, start = sp_Cpos, end = -1)
    
    out_obj <- SignalpResult(in_fasta = fasta,
                             out_fasta = out_fasta_sp, 
                             mature_fasta = cropped_fasta, 
                             sp_version = version,
                             sp_tibble = sp)
    if (validObject(out_obj)) {return(out_obj)}
    
  } else {
    message('Allowed versions include...:')
    message(cat(allowed_versions))
    message('Allowed organisms iclude...')
    message(cat(allowed_organisms))}
    stop('Input signalp version or specified organism type are invalid.')
}


# ### generate out_fasta
r3 <- signalp(w, version = 3, 'euk')
r4 <- signalp(w, version = 4, 'euk')
r2 <- signalp(w, version = 2, 'euk') # does not return correct outout, repeted ids
