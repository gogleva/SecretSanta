#' parse_signalp function
#'
#' This function parses signalp2 and signalp3 output
#' @param input signalp2 or 3 input
#' @param input_type  path: file with text output from signalp2 or signalp3 \cr 
#'                    system_call: direct output from signalp2/3 system call
#' @export
#' @examples
#' res2_path <- parse_signalp(input = "/home/anna/anna/Labjournal/SecretSanta/inst/extdata/sample_prot_signalp2_out", input_type = "path")
#' con <- system("/home/anna/anna/Labjournal/SecretSanta_external/signalp-2.0/signalp -t euk SecretSanta/inst/extdata/sample_prot.fasta", intern = TRUE)
#' res2_system <- parse_signalp(input = con, input_type = "system_call")

parse_signalp <- function(input, input_type) {
  # helper function for gene ids
  clean_geneids <- function(x) {gsub('>', '', unlist(stringr::str_split(x, " "))[1])}
  # helper function for C-score, Y-score and S-score: split line with varibale number of spaces
  clean_score <- function(x) {as.numeric(stringr::str_split(x, "\\s+")[[1]][c(4,5 )])}
  # helper function for S mean
  clean_mean <- function(x) {strsplit(x, "\\s+")[[1]][c(4,5)]}
  # helper fucntion for prediction result:
  clean_status <- function(x) {gsub('Prediction: ', '', x)}
  # read data
  if (input_type == 'path') {
    data <- readLines(input)
  } else if (input_type == 'system_call'){
    data <- input #system call already captured in a character object
  }
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
  res <- tibble::as.tibble(data.frame(gene_ids_fixed,
                              t(max_C_fixed),
                              t(max_Y_fixed),
                              t(max_S_fixed),
                              t(mean_S_fixed),
                              Status_fixed))
  names(res) <- c("gene_id", "Cpos", "Cmax",
                  "Ypos", "Ymax", "Spos",
                  "Smax", "Srange", "Smean", "Prediction")
  #filter entries predicted to contain signal peptide
  return(res %>% filter(Prediction == 'Signal peptide'))
}

###tests

#signalp2 output
#res2 <- parse_signalp("/home/anna/anna/Labjournal/SecretSanta/inst/extdata/sample_prot_signalp2_out")
#signalp3 output
#res3 <- parse_signalp("SecretSanta/inst/extdata/sample_prot_signalp3_out")


# # problem is in the parserL signalp3 output is ok
con3 <- system("/home/anna/anna/Labjournal/SecretSanta_external/signalp-3.0/signalp -t euk SecretSanta/inst/extdata/sample_prot.fasta", intern = TRUE)
res3_system <- parse_signalp(input = con3, input_type = "system_call")

# # signalp2 output is tropbled - fix this tomorrow
con2 <- system("/home/anna/anna/Labjournal/SecretSanta_external/signalp-2.0/signalp -t euk SecretSanta/inst/extdata/sample_prot.fasta", intern = TRUE)
res2_system <- parse_signalp(input = con2, input_type = "system_call")

res2_path <- parse_signalp(input = "/home/anna/anna/Labjournal/SecretSanta/inst/extdata/sample_prot_signalp2_out", input_type = "path")

#identical(con2, con3)
# => outputs of signalp2 and signalp3 are different, sp2 parsing fails



