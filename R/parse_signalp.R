#' parse_signalp function
#'
#' This function parses signalp2 and signalp3 output
#' @param signalp_output  text output from signalp2 or signalp3, test: read from the file
#' @export
#' @examples
#' parse_signalp(some_param)

parse_signalp <- function(signalp_output) {
  clean_geneids <- function(x) {gsub('>', '', unlist(str_split(x, " "))[1])}
  data <-  readLines("SecretSanta/inst/extdata/sample_prot_signalp2_out")
  gene_ids <- data[(grep("# Measure  Position  Value  Cutoff  signal peptide?", data) - 1)]
  gene_ids_fixed <- (sapply(gene_ids, gene_id_rep, USE.NAMES = FALSE))
  
}


data <- readLines("SecretSanta/inst/extdata/sample_prot_signalp2_out") #now read from the file, later - update this
# raw lines

max_C <- data[grep("max. C", data)]  
max_Y <- data[grep("max. Y", data)]
max_S <- data[grep("max. S", data)]
mean_S <- data[grep("mean S", data)]

#clean them!



