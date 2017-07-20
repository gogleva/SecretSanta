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
  # helper function for C-score
  clean_maxC <- function(x) {as.numeric(unlist(str_split(x, "  "))[c(3,6)])}
  # read data
  data <-  readLines("SecretSanta/inst/extdata/sample_prot_signalp2_out")
  # extract gene ids
  gene_ids <- data[(grep("# Measure  Position  Value  Cutoff  signal peptide?", data) - 1)]
  gene_ids_fixed <- (sapply(gene_ids, clean_geneids, USE.NAMES = FALSE))
  # extract max C score and position
  max_C <- data[grep("max. C", data)]  
  max_C_fixed <- (sapply(max_C, clean_maxC, USE.NAMES = FALSE))
}


data <- readLines("SecretSanta/inst/extdata/sample_prot_signalp2_out") #now read from the file, later - update this
# raw lines

max_C <- data[grep("max. C", data)]  
max_Y <- data[grep("max. Y", data)]
max_S <- data[grep("max. S", data)]
mean_S <- data[grep("mean S", data)]

#clean them!



