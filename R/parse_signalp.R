#' convert output of Signalp-2.0 and Signalp-3.0 to Signalp-4.0++ format
#'
#' This function parses the output of the command line tools \code{signalp2} and \code{signalp3} to standardize outputs for data analysis.\cr
#' \cr
#' Alternatively, \code{parse_signalp} can be called independently on outputs of
#' \code{signalp2} and \code{signalp3} and captured in a system call or stored in a file.
#' @param input output of the command line tools \code{signalp2} or \code{signalp3}/
#' @param input_type a character string with the following options:
#' \itemize{
#' \item \code{input_type = "path"}  - path to a file with text output from \code{signalp2}
#' or \code{signalp3}
#' \item \code{input_type = "system_call"} - output from the \code{signalp2}
#' or \code{signalp3} system call
#' }
#' @param method which prediction method to use. Options are:
#' \itemize{
#' \item \code{method = "hmm"} - for HMM-based predictions
#' \item \code{method = "nn"} - for prediction based on neural networks
#' }
#' @param pred_filter filter for the type of prediction. Options are:
#' \itemize{
#' \item \code{pred_filter = "Signal peptide"}
#' \item \code{pred_filter = "Signal anchor"}
#' \item \code{pred_filter = "Non-secretory protein"}
#' \item \code{pred_filter = "all"} - in case all three filter options shall be included
#' }
#' @param version version of SignalP used: 2.0 or 3.0
#' @param source_fasta source fasta file, required to resque names when signalp2, nn method is used
#' @return parsed \code{signalp2}
#' or \code{signalp3} output, organised in a \code{\link[tibble]{tibble}} object.
#' @export
#' @examples
#' # Example 1: parse signalp2 output, stored in a file:
#' s_path <- system.file("extdata", "sample_prot_sp2_hmm_out",
#' package = "SecretSanta")
#' parse_signalp(input = s_path, input_type = "path", pred_filter = "Signal peptide", version = 3, method = 'hmm')
#'
#' # alternatively users can select for all prediction filters
#' parse_signalp(input = s_path, input_type = "path", pred_filter = "all", method = 'hmm', version = 2)
#'
#' # Example2: parse signalp2 output
#' # captured in a system call. Note, here we assume that
#' # signalp2 is accessible via $PATH.
#' s_fasta <- system.file("extdata", "small_prot.fasta",
#' package = "SecretSanta")
#' # capture system call:
#' con_hmm <- system(paste('signalp2 -t euk -f short -m hmm -trunc 70', s_fasta), intern = TRUE)
#' con_nn <- system(paste('signalp3 -t euk -f short -m nn -trunc 70', s_fasta), intern = TRUE)
#' # parse captured system call:
#' parse_signalp(input = con_hmm, input_type = "system_call", method = 'hmm', version = 2)
#' parse_signalp(input = con_nn, input_type = "system_call",
#' method = 'nn', version = 3, pred_filter = 'all')
#' @seealso \code{\link{signalp}}

parse_signalp <-
    function(input,
             input_type = c('path', 'system_call'),
             method = c('nn', 'hmm'),
             pred_filter = "Signal peptide",
             version = c(2, 3),
             source_fasta = NULL
             ){

        input_type <- match.arg(input_type)
        
        if (missing(method)) {
            stop('missing argument: method')
        }
        method <- match.arg(method)
        
        if (missing(version)) {
            stop('missing argument: version')
        }

        # ----- source_fasta set to default or complain if not provided
        # when sp2+nn
        
        if (all(version == 2, method == 'nn')) {
            if (is.null(source_fasta)) {
                stop("please provide source_fasta")
            } else {
                if (!(file.exists(source_fasta))) {
                    stop('source_fasta file does not exist, please check the supplied file name')
                }
            }
        }

        if (!is.element(input_type, c('path', 'system_call')))
          stop("please specify either input_type = 'path' or input_type = 'system_call'.", call. = FALSE)

        if (!is.element(pred_filter, c("Signal peptide", "Signal anchor", "Non-secretory protein", "all")))
          stop("please provide a valid filter type ...", call. = FALSE)

        message("signalp output is imported and filter '", pred_filter,"' is applied ...")
        
        # read data
        if (input_type == 'path') {
            data <- read.table(input)
        } else if (input_type == 'system_call') {
            tmp_con <- tempfile()
            write(input, tmp_con) # put everything in a tmp file
            data <- read.table(tmp_con)
        }
        
        # testing:
        inp2_nn <- system(paste('signalp2 -t euk -f short -m nn -trunc 70', s_fasta), intern = TRUE)
        inp3_nn <- system(paste('signalp3 -t euk -f short -m nn -trunc 70', s_fasta), intern = TRUE)   
        
        # proceed differently, depending on a prediction method selected
        if (method == 'nn') {
            
            if (version == 2) {
                # we need to rescue seqeunce names, nn method crops
                # them too much and can create duplicate gene_ids when
                # siganlp2 is used
                
                names(data) <- c('gene_id', 'Cmax', 'Cpos',
                                 'C_pred', 'Ymax', 'Ypos', 'Y_pred',
                                 'Smax', 'Spos', 'S_pred', 'Smean',
                                 'S_predm')
                
                data$Prediction_YN <- data$Ymax 
                
                fasta <- readAAStringSet(source_fasta)
                cropped_names <- unname(sapply(names(fasta), crop_names))
                data$gene_id <- cropped_names
                
            } else {
                
           names(data) <- c('gene_id', 'Cmax', 'Cpos',
                             'C_pred', 'Ymax', 'Ypos', 'Y_pred',
                             'Smax', 'Spos', 'S_pred', 'Smean',
                             'S_predm', 'D', 'Prediction_YN'
                             )
            }
            
            data$Sprob <- data$Prediction <- NA
            data <- as.tibble(data) %>% 
                dplyr::mutate(
                    Prediction = ifelse(Prediction_YN == 'Y',
                                        "Signal peptide", NA))
            
            
        } else if (method == 'hmm') {
            # narrow table
            names(data) <- c('gene_id', 'Prediction', 'Cmax',
                             'Cpos', 'v5', 'Sprob', 'Prediction_YN')
            # columns absent from the hmm-output
            data$Ymax <- data$Ypos <- data$Smax <-
                data$Spos <- data$Smean <- NA
            
            # replace SQA codes to be consistent with nn output
            data <- as_tibble(data) %>% 
                dplyr::mutate(Prediction = case_when(
                Prediction == 'S' ~ "Signal peptide",
                Prediction == 'A' ~ "Signal anchor",
                Prediction == 'Q' ~ "Non-secretory protein"))
        }
        
        # check that there are no duplicated gene ids:
        if (any(duplicated(data$gene_id))) {
            stop('gene_ids vector contains duplicated elements ...', 
                 call. = FALSE)
        }
        
        # reformat table to be compatible with siganlp4++ output
        
        gene_id <- Cmax <- Cpos <- Ymax <- Ypos <- Smax <- NULL
        Spos <- Smean <- Prediction <- Prediction_YN <- NULL
        res <- dplyr::select(data, gene_id, Cmax,
                             Cpos, Ymax, Ypos, Smax,
                             Spos, Smean, Prediction)             
        #filter entries predicted to contain signal peptide
        if (pred_filter != "all")
        res <- dplyr::filter(res, Prediction == pred_filter)
        if (pred_filter == "all")
        res <- dplyr::filter(res, Prediction %in% c("Signal peptide", "Signal anchor", "Non-secretory protein"))

        #Smean to numeric
        res$Smean <- as.numeric(as.character(res$Smean))
        message("import completed!")
        return(res)
    }