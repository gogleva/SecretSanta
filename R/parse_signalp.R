#' parse_signalp function
#'
#' This function parses signalp2 and signalp3 output
#' @param signalp_output  text output from signalp2 or signalp3, test: read from the file
#' @export
#' @examples
#' parse_signalp(some_param)

parse_signalp <- function(signalp_output) {
  # helper function for gene ids
  clean_geneids <- function(x) {gsub('>', '', unlist(str_split(x, " "))[1])}
  # helper function for C-score, Y-score and S-score: split line with varibale number of spaces
  clean_score <- function(x) {as.numeric(strsplit(x, "\\s+")[[1]][c(4,5)])}
  # read data
  data <-  readLines("SecretSanta/inst/extdata/sample_prot_signalp2_out")
  # extract gene ids
  gene_ids <- data[(grep("# Measure  Position  Value  Cutoff  signal peptide?", data) - 1)]
  gene_ids_fixed <- (sapply(gene_ids, clean_geneids, USE.NAMES = FALSE)) #fin
  # extract max C score and position
  max_C_fixed <- sapply(data[grep("max. C", data)], clean_score, USE.NAMES = FALSE) #fin
  # extract max Y score and position
  max_Y_fixed <- sapply(data[grep("max. Y", data)], clean_score, USE.NAMES = FALSE) #fin
  # extract max S score and position
  max_S_fixed <- sapply(data[grep("max. S", data)], clean_score, USE.NAMES = FALSE) #fin
  # extract mean S score and position
#  mean_S_fixed <- sapply(data[grep("mean S", data)], clean_score, USE.NAMES = FALSE) #fin
}
  
 
#split string with variable number of spaces:
strsplit(max_C[10], "\\s+")[[1]]

data <- readLines("SecretSanta/inst/extdata/sample_prot_signalp2_out") #now read from the file, later - update this
# raw lines

max_Y <- data[grep("max. Y", data)]
max_S <- data[grep("max. S", data)]
mean_S <- data[grep("mean S", data)]

