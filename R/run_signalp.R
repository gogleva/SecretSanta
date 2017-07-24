#' signalp function
#'
#' This function calls local SignalP
#' Please ensure that respective version of SignalP is downloaded, installed and respective path is added to $PATH variable.
#' For details see ... explain with aliases and examples
#' @param proteins input file with proteins
#' @param version  specify SignalP version to run; versions available: \cr
#'                 2.0 - http://www.cbs.dtu.dk/services/SignalP-2.0/ \cr
#'                 3.0 - http://www.cbs.dtu.dk/services/SignalP-3.0/ \cr
#'                 4.1 - http://www.cbs.dtu.dk/services/SignalP-4.1/ \cr
#' @param organism_type euk, gram+, gram-                 
#' @param output_type options: tibble/table with results; - this is done \cr
#'                    fasta file with candidate secreted peptides, full length - not done; \cr
#'                    fasta file with mature sequences - not done
#' @export
#' @examples
#' result <- signalp(proteins = "SecretSanta/inst/extdata/sample_prot.fasta", organism_type = 'euk', version = 4)

# alterantively: capture all possible outputs in a complex object, like DEseq objects

signalp <- function(proteins, version, organism_type, output_type) {
  message("running signalP locally...")
  # check input fasta files for illegal characters and other incompatibilities
  
  fasta <- readLines(proteins)
  if (any(grepl('[*$]',fasta[!grepl('>', fasta)]))){
    stop('input fasta file cotains stop codons "*", please remove them and re-run')
  }else{
    # determine which of signalp paths we need to use:
    signalp_version <- paste("signalp", version, sep = '')
    message(signalp_version)
    full_pa <- as.character(secret_paths %>% filter(tool == signalp_version) %>% select(path))
  
    if (version >= 4) {
      # running signalp versions 4 and 4.1
      result <- tibble::as.tibble(read.table(text = (system(paste(full_pa, "-t", organism_type, proteins), intern = TRUE))))
      names(result) <- c("gene_id", "Cmax", "Cpos",
                       "Ymax", "Ypos", "Smax",
                       "Spos", "Smean", "D",
                       "Prediction", "Dmaxcut", "Networks-used")
      # returns tibble for candidate secreted proteins only
      return(result %>% filter(Prediction == 'Y'))
    } else if(version < 4) {
      # running signalp versions 2 and 3, call parser for the initial output
      message('ancient signalp, calling parser for the output...')
      con <- system(paste(full_pa, "-t", organism_type, proteins), intern = TRUE)
      result <- parse_signalp(input = con, input_type = "system_call")
    }
  }
}
