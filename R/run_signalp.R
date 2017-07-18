#' signalp function
#'
#' This function calls local SignalP
#' Please ensure that respective version of SignalP is downloaded, installed and respective path is added to $PATH variable.
#' For details see ... explain with aliases and examples
#' @param proteins input file/data frame with proteins
#' @param version  specify SignalP version to run; versions available: \cr
#'                 2.0 - http://www.cbs.dtu.dk/services/SignalP-2.0/ \cr
#'                 3.0 - http://www.cbs.dtu.dk/services/SignalP-3.0/ \cr
#'                 4.1 - http://www.cbs.dtu.dk/services/SignalP-4.1/ \cr
#' @param organism_type euk, gram+, gram-                 
#' @export
#' @examples
#' result <- signalp(proteins = "SecretSanta/inst/extdata/sample_prot.fasta", organism_type = 'euk', version = 4)


signalp <- function(proteins, version, organism_type) {
  message("running signalP locally...")
  #determine which of signalp paths we need to use:
  signalp_version <- paste("signalp", version, sep = '')
  print(signalp_version)
  full_pa <- as.character(secret_paths %>% filter(tool == signalp_version) %>% select(path))
  if (version >= 4) {
    #running signalp versions 4 and 4.1
    result <- tibble::as.tibble(read.table(text = (system(paste(full_pa, "-t", organism_type, proteins), intern = TRUE))))
    names(result) <- c("gene_id", "Cmax", "Cpos",
                       "Ymax", "Ypos", "Smax",
                       "Spos", "Smean", "D",
                       "Status", "Dmaxcut", "Networks-used")
    # returns tibble for candidate secreted proteins only
    return(result %>% filter(Status == 'Y'))
  } else {
    print('ancient signalp, requires parser first')
  }
}

