#' parse_signalp function
#' 
#' This function parses signalp2 and signalp3 output and is called internally
#' in the \code{\link{signalp}} function to standardize outputs.\cr
#' \cr
#' Alternatively, parse_signalp can be called independently on outputs of 
#' signalp2 and signalp3 captured in a system call or stored in a file.
#' @param input output of signalp2 or signalp3
#' @param input_type
#' \strong{path}  path to a file with text output from signalp2 or signalp3;\cr
#' \strong{system_call} output from signalp2/3 system call;
#' @return parsed signalp2/3 output, organised in a tibble object.
#' @export
#' @examples
#' # Example 1: parse signalp2 output, stored in a file:
#' s_path <- system.file("extdata",
#'                       "sample_prot_signalp2_out2",
#'                        package = "SecretSanta") 
#' parse_signalp(input = s_path,
#'               input_type = "path")
#' 
#' # Example2: parse signalp2 output
#' # captured in a system call. Note, here we assume that
#' # signalp2 is accessible via $PATH.
#' s_fasta <- system.file("extdata",
#'                        "small_prot.fasta",
#'                         package = "SecretSanta") 
#' # capture system call:
#' con <- system(paste('signalp2 -t euk', s_fasta),
#'               intern = TRUE)
#' 
#' # parse captured system call:
#' parse_signalp(input = con,
#'               input_type = "system_call")

parse_signalp <-
  function(input,
           input_type = c('path', 'system_call')) {
    
    input_type <- match.arg(input_type)
    
    # helper function to rescue gene ids
    clean_geneids <- function(x) {
      gsub('>', '', unlist(stringr::str_split(x, " "))[1])
    }
    
    # helper function to parse C-score, Y-score and S-score:
    # split line with varibale number of spaces
    clean_score <- function(x) {
      as.numeric(stringr::str_split(x, "\\s+")[[1]][c(4, 5)])
    }
    
    # helper function to get signal peptide cleavage site (from HMM prediction)
    # NN predictions often output wrong coordinates
    clean_cleavege <- function(x) {
      as.numeric(tail(unlist(stringr::str_split(x, "\\s+")), n = 1))
    }
    
    # helper function to parse S mean
    clean_mean <- function(x) {
      strsplit(x, "\\s+")[[1]][c(4, 5)]
    }
    
    # helper function to parse prediction summary:
    clean_status <- function(x) {
      gsub('Prediction: ', '', x)
    }
    
    # read data
    if (input_type == 'path') {
      data <- readLines(input)
    } else if (input_type == 'system_call') {
      data <- input #system call already captured in a character object
    }
    
    # extract gene ids
    gene_ids <- data[(grep("SignalP-HMM result:", data) + 1)]
    gene_ids_fixed <-
      (sapply(gene_ids, clean_geneids, USE.NAMES = FALSE))
    
    # check that there are no duplicated gene ids:
    if (any(duplicated(gene_ids))) {
      stop('gene_ids vector contains duplicated elements')
    }
    
    # extract max C score and position
    max_C_fixed <- sapply(data[grep("max. C", data)],
                          clean_score,
                          USE.NAMES = FALSE)
    
    # extract max Y score and position
    C_pos <-
      sapply(data[grep("Max cleavage site probability:", data)],
             clean_cleavege,
             USE.NAMES = FALSE)
    
    max_Y_fixed <- sapply(data[grep("max. Y", data)],
                          clean_score,
                          USE.NAMES = FALSE)
    
    # extract max S score and position
    max_S_fixed <- sapply(data[grep("max. S", data)],
                          clean_score,
                          USE.NAMES = FALSE)
    
    # extract mean S score and position
    mean_S_fixed <- sapply(data[grep("mean S", data)],
                           clean_mean,
                           USE.NAMES = FALSE)
    
    # extract prediction summary
    Status_fixed <- sapply(data[grep("Prediction: ", data)],
                           clean_status,
                           USE.NAMES = FALSE)
    
    # assemble result object
    res <- tibble::as.tibble(data.frame(
      gene_ids_fixed,
      t(max_C_fixed)[, 2],
      C_pos,
      t(max_Y_fixed),
      t(max_S_fixed),
      t(mean_S_fixed),
      Status_fixed
    ))
    
    names(res) <- c(
      "gene_id",
      "Cmax",
      "Cpos",
      "Ypos",
      "Ymax",
      "Spos",
      "Smax",
      "Srange",
      "Smean",
      "Prediction"
    )
    
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
    res <- res %>% filter(res$Prediction == 'Signal peptide')
    
    #Smean to numeric
    res$Smean <- as.numeric(as.character(res$Smean))
    
    return(res)
  }
