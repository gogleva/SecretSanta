#' parse_signalp function
#' 
#' This function parses signalp2 and signalp3 output. Is used internally in signalp function; can be called independently on outputs of signalp2/3 tools captured in a system call or stored in a file.
#' @param input full output of signalp2 or signalp3 call
#' @param input_type  specify input type: 'path' or 'system_call'
#' \itemize{
#' \item 'path' - path to a file with text output from signalp2 or signalp3;
#' \item 'system_call' - direct output from signalp2/3 system call;
#' }
#' @return parsed signalp2/3 output, organised in a tibble object.
#' @export
#' @examples
#' 
#' # Parse signalp2 output, stored in a file:
#' s_path <- system.file("extdata", "sample_prot_signalp2_out", package = "SecretSanta") 
#' parse_sp_path <- parse_signalp(input = s_path, input_type = "path")
#' 
#' # Parse signalp2 output, obtained from a system call:
#' s_fasta <- system.file("extdata", "sample_prot.fasta", package = "SecretSanta") 
#' secret_paths <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
#' sp2_path <- secret_paths %>% filter(tool == 'signalp2') %>% select(path)
#' 
#' # capture system call:
#' con <- system(paste(sp2_path, '-t euk', s_fasta), intern = TRUE)
#' 
#' # parse captured system call:
#' parse_signalp(input = con, input_type = "system_call")

parse_signalp <- function(input, input_type) {
  # helper function for gene ids
  clean_geneids <- function(x) {gsub('>', '', unlist(stringr::str_split(x, " "))[1])}
  # helper function for C-score, Y-score and S-score: split line with varibale number of spaces
  clean_score <- function(x) {as.numeric(stringr::str_split(x, "\\s+")[[1]][c(4,5 )])}
  # helper function to get signal peptide cleavage site
  clean_cleavege <- function(x) {as.numeric(tail(unlist(stringr::str_split(x, "\\s+")), n = 1))}
  # helper function for S mean
  clean_mean <- function(x) {strsplit(x, "\\s+")[[1]][c(4,5)]}
  # helper function for prediction summary:
  clean_status <- function(x) {gsub('Prediction: ', '', x)}
  
  # read data
  if (input_type == 'path') {
    data <- readLines(input)
  } else if (input_type == 'system_call'){
    data <- input #system call already captured in a character object
  }
  # extract gene ids
  gene_ids <- data[(grep("SignalP-HMM result:", data) + 1)]
  gene_ids_fixed <- (sapply(gene_ids, clean_geneids, USE.NAMES = FALSE))
  
  # check that there are no duplicated gene ids:
  if (identical(length(gene_ids), length(unique(gene_ids)))){} else {stop('gene_ids vector contains duplicated elements')}
  
  # extract max C score and position
  max_C_fixed <- sapply(data[grep("max. C", data)], clean_score, USE.NAMES = FALSE)
  # extract max Y score and position
  C_pos <- sapply(data[grep("Max cleavage site probability:", data)], clean_cleavege, USE.NAMES = FALSE)
  
  
  max_Y_fixed <- sapply(data[grep("max. Y", data)], clean_score, USE.NAMES = FALSE)
  # extract max S score and position
  max_S_fixed <- sapply(data[grep("max. S", data)], clean_score, USE.NAMES = FALSE)
  # extract mean S score and position
  mean_S_fixed <- sapply(data[grep("mean S", data)], clean_mean, USE.NAMES = FALSE)
  Status_fixed <- sapply(data[grep("Prediction: ", data)], clean_status, USE.NAMES = FALSE)
  res <- tibble::as.tibble(data.frame(gene_ids_fixed,
                              t(max_C_fixed)[2],
                              C_pos,
                              t(max_Y_fixed),
                              t(max_S_fixed),
                              t(mean_S_fixed),
                              Status_fixed))
  names(res) <- c("gene_id", "Cmax", "Cpos",
                  "Ypos", "Ymax", "Spos",
                  "Smax", "Srange", "Smean", "Prediction")

  #re-order columns to match signalp4 output
  
  res <- res %>% select("gene_id", 
                        "Cmax",
                        "Cpos",
                        "Ymax",
                        "Ypos",
                        "Smax",
                        "Spos",
                        "Smean",
                        "Prediction")
  
  #filter entries predicted to contain signal peptide
  return(res %>% filter(Prediction == 'Signal peptide'))
}

