#' convert output of Signalp-2.0 and Signalp-3.0 to Signalp-4.0++ format
#' 
#' This function parses signalp2 and signalp3 output and is called internally
#' in the \code{\link{signalp}} function to standardize outputs.\cr
#' \cr
#' Alternatively, parse_signalp can be called independently on outputs of 
#' signalp2 and signalp3 captured in a system call or stored in a file.
#' @param input output of signalp2 or signalp3
#' @param input_type
#' \strong{path}  path to a file with text output from signalp2 
#' or signalp3;\cr
#' \strong{system_call} output from signalp2/3 system call;
#' @return parsed signalp2/3 output, organised in a tibble object.
#' @export
#' @examples
#' # Example 1: parse signalp2 output, stored in a file:
#' s_path <- system.file("extdata", "sample_prot_signalp2_out2",
#' package = "SecretSanta") 
#' parse_signalp(input = s_path, input_type = "path")
#' 
#' # Example2: parse signalp2 output
#' # captured in a system call. Note, here we assume that
#' # signalp2 is accessible via $PATH.
#' s_fasta <- system.file("extdata", "small_prot.fasta", 
#' package = "SecretSanta") 
#' # capture system call:
#' con <- system(paste('signalp2 -t euk', s_fasta), intern = TRUE)
#' 
#' # parse captured system call:
#' parse_signalp(input = con, input_type = "system_call")

parse_signalp <-
    function(input,
             input_type = c('path', 'system_call')) {
        input_type <- match.arg(input_type)
        
        # helper functions----
        # helper function to rescue gene ids
        clean_geneids <- function(x) {
            gsub('>', '', unlist(strsplit(x, " "))[1])
        }
       
        # helper function to parse C-score, Y-score and S-score:
        # split line with varibale number of spaces
        clean_score <- function(x) {
            as.numeric(strsplit(x, "\\s+")[[1]][c(4, 5)])
        }
        
        # helper function to get signal peptide cleavage site (from HMM prediction)
        # NN predictions often output wrong coordinates
        clean_cleavege <- function(x) {
            as.numeric(tail(unlist(strsplit(x, "\\s+")), n = 1))
        }
        
        # helper function to parse S mean
        clean_mean <- function(x) { strsplit(x, "\\s+")[[1]][c(4, 5)]}
        
        # helper function to parse prediction summary:
        clean_status <- function(x) { gsub('Prediction: ', '', x) }
        
        # end of helper functions----
        
        # read data
        if (input_type == 'path') {
            data <- readLines(input)
        } else if (input_type == 'system_call') {
            data <- input #system call already captured in a character object
        }
        
        # extract gene ids
        gene_ids <- data[(grep("SignalP-HMM result:", data) + 1)]
        gene_ids_fixed <- (sapply(gene_ids, clean_geneids, USE.NAMES = FALSE))
        
        # check that there are no duplicated gene ids:
        if (any(duplicated(gene_ids))) {
            stop('gene_ids vector contains duplicated elements')
        }
        
        # clean the remaining fields in a more optimal/functional way
        
        
        # TESTING:: data <- system(paste('signalp2 -t euk', s_fasta), intern = TRUE)
        
        # organise cleaning functions in a list:
        clean_stats_fun <- list(clean_cleavege, clean_mean, clean_status)
            
#            list(clean_score, clean_cleavege, clean_mean,
 #                        clean_status)
        
        # helper function to grep plain text data by mathcing patterns
        grep_data <- function(grep_expr) {data[grep(grep_expr, data)]}
        
        # key hooks to grep the lines with required info
        grep_param <- c("max. C", "max. Y", "max. S",
                        "Max cleavage site probability:", "mean S",
                        "Prediction: ")
        
        arg_list <- lapply(grep_param, grep_data)
        
        # 2 separte maps - one to map a list of functions over a list of lines;
        # the other - just to map clean_score over a list of lines and transpose
        # the output
        
        map_stats <- Map(function(x,y) sapply(y, x, USE.NAMES = FALSE),
                         clean_stats_fun, arg_list[4:6])
        names(map_stats) <- c("C_pos", "mean_S_fixed", "Status_fixed")
        
        map_scores <- Map(function(x) t(sapply(x, clean_score, USE.NAMES = FALSE)),
                          arg_list[1:3])    
        names(map_scores) <- c("C_max_fixed", "Y_max_fixed", "S_max_fixed")
        
       # initial combined map function
       # mymap <- Map(function(x,y) sapply(y,x, USE.NAMES = FALSE),
        #             fun_fun_list, arg_list)
        

        ## This is not super elegant, but works - apply transposition better
        res <- tibble::as_tibble(data.frame(
            gene_ids_fixed,
            t(mymap[[1]])[,2],
            mymap[[2]],
            t(mymap[[3]]),
            t(mymap[[4]]),
            t(mymap[[5]]),
            mymap[[6]]
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
        # MAY BE I COULD AVOID THIS STEP AT ALL?
        res <- res[c("gene_id", "Cmax", "Cpos",
                     "Ymax", "Ypos", "Smax",
                     "Spos", "Smean", "Prediction")]
        
        
        #filter entries predicted to contain signal peptide
        res <- res[res$Prediction == 'Signal peptide',]
        
        #Smean to numeric
        res$Smean <- as.numeric(as.character(res$Smean))
        return(res)
    }
