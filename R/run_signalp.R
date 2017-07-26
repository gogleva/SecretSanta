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
#' r <- signalp(w, version = 4, 'euk')
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
    out <- tempfile() #create a temporary file for fasta
    Biostrings::writeXStringSet(fasta, out) #write tmp fasta file
    # make a system call of signalp based on the tmp file
    signalp_version <- paste("signalp", version, sep = '')
    message(signalp_version)
    full_pa <- as.character(secret_paths %>% filter(tool == signalp_version) %>% select(path))
    result <- tibble::as.tibble(read.table(text = (system(paste(full_pa, "-t", organism_type, out), intern = TRUE))))
    return(result)
  } else {
    message('Allowed versions include...:')
    message(cat(allowed_versions))
    message('Allowed organisms iclude...')
    message(cat(allowed_organisms))}
    stop('Input signalp version or specified organism type are invalid.')
}

  #   if (version >= 4) {
  #     # running signalp versions 4 and 4.1
  #     result <- tibble::as.tibble(read.table(text = (system(paste(full_pa, "-t", organism_type, proteins), intern = TRUE))))
  #     names(result) <- c("gene_id", "Cmax", "Cpos",
  #                      "Ymax", "Ypos", "Smax",
  #                      "Spos", "Smean", "D",
  #                      "Prediction", "Dmaxcut", "Networks-used")
  #     # returns tibble for candidate secreted proteins only
  #     return(result %>% filter(Prediction == 'Y'))
  #   } else if(version < 4) {
  #     # running signalp versions 2 and 3, call parser for the initial output
  #     message('ancient signalp, calling parser for the output...')
  #     con <- system(paste(full_pa, "-t", organism_type, proteins), intern = TRUE)
  #     result <- parse_signalp(input = con, input_type = "system_call")
  #   }
  # }
}


## test inputs for system calls

bb
# test system call with file input
f_input <- "/home/anna/anna/Labjournal/SecretSanta/inst/extdata/sample_prot.fasta"
s_con <- system(paste("/home/anna/anna/Labjournal/SecretSanta_external/signalp-4.1/signalp -t euk", f_input))


# what we really want - to use and object as an input

bb_seqs <- as.character(bb, use.names = TRUE)
bb_seqs2 <-toString(bb, use.names = TRUE)

ss_tmp <- "/home/anna/anna/Labjournal/SecretSanta_tmp/"
writeXStringSet(bb, file = paste(ss_tmp, "bb_tmp.fasta", sep = ''), format = 'fasta', width=80)

out1 <- tempfile()
writeXStringSet(bb, out1)

#try system call for temporary file
filepath <- system.file("extdata", "someORF.fa", package="Biostrings")
fasta.info(filepath)
x <- readDNAStringSet(filepath)
out1 <- tempfile()
s_con_tmp <- system(paste("/home/anna/anna/Labjournal/SecretSanta_external/signalp-4.1/signalp -t euk", out1))
#=> it works, now need to implement it in the signalp function

#
tt <- 'test branching with package'

