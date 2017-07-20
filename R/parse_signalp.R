#' parse_signalp function
#'
#' This function parses signalp2 and signalp3 output
#' @param data_path  text output from signalp2 or signalp3, test: read from the file
#' @export
#' @examples
#' parse_signalp("SecretSanta/inst/extdata/sample_prot_signalp2_out")

parse_signalp <- function(data_path) {
  # helper function for gene ids
  clean_geneids <- function(x) {gsub('>', '', unlist(str_split(x, " "))[1])}
  # helper function for C-score, Y-score and S-score: split line with varibale number of spaces
  clean_score <- function(x) {as.numeric(strsplit(x, "\\s+")[[1]][c(4,5 )])}
  # helper function for S mean
  clean_mean <- function(x) {strsplit(x, "\\s+")[[1]][c(4,5)]}
  # helper fucntion for prediction result:
  clean_status <- function(x) {gsub('Prediction: ', '', x)}
  # read data
  data <-  readLines(data_path)
  # extract gene ids
  gene_ids <- data[(grep("# Measure  Position  Value  Cutoff  signal peptide?", data) - 1)]
  gene_ids_fixed <- (sapply(gene_ids, clean_geneids, USE.NAMES = FALSE))
  # extract max C score and position
  max_C_fixed <- sapply(data[grep("max. C", data)], clean_score, USE.NAMES = FALSE)
  # extract max Y score and position
  max_Y_fixed <- sapply(data[grep("max. Y", data)], clean_score, USE.NAMES = FALSE)
  # extract max S score and position
  max_S_fixed <- sapply(data[grep("max. S", data)], clean_score, USE.NAMES = FALSE)
  # extract mean S score and position
  mean_S_fixed <- sapply(data[grep("mean S", data)], clean_mean, USE.NAMES = FALSE)
  Status_fixed <- sapply(data[grep("Prediction: ", data)], clean_status, USE.NAMES = FALSE)
  res <- as.tibble(data.frame(gene_ids_fixed,
                              t(max_C_fixed),
                              t(max_Y_fixed),
                              t(max_S_fixed),
                              t(mean_S_fixed),
                              Status_fixed))
  names(res) <- c("gene_id", "Cpos", "Cmax",
                  "Ypos", "Ymax", "Spos",
                  "Smax", "Srange", "Smean", "Status")
  #filter entries predicted to contain signal peptide
  return(res %>% filter(Status == 'Signal peptide'))
}
  
res2 <- parse_signalp("SecretSanta/inst/extdata/sample_prot_signalp2_out")

data <- readLines("SecretSanta/inst/extdata/sample_prot_signalp2_out")
