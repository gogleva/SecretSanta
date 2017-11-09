#' crop_names function
#'
#' helper function crop long names for AAStringSet object, return
#' character vector
#' @return character vector 
#' @keywords internal
## @export 

crop_names <- function(x){unlist(stringr::str_split(x, " "))[1]}

#' truncate_seq
#' 
#' helper function to truncate long sequences or throw them away, otherwise
#' signalp will break (at least signalp2 and signalp3 will)
#' @return truncated AAStringSet
#' @param truncate truncate or not
#' @param seq_set XStringSet
#' @param threshold length threshold
#' @keywords internal
## @export 

truncate_seq <- function(truncate, seq_set, threshold) {
    drop_n <- length(seq_set[width(seq_set) >= threshold])

    if (drop_n == 0)
        return(seq_set)

    if (truncate == FALSE) {
        seq_set <- seq_set[width(seq_set) < threshold]
        warning(paste(drop_n, 'long sequenses have been thrown away'))
        return(seq_set)

    } else if (truncate == TRUE) {
        message(paste(drop_n, 'sequences to be truncated'))
        seq_keep <-
            seq_set[width(seq_set) < threshold] # not so long sequences
        seq_trunc <-
            seq_set[width(seq_set) >= threshold] # sequences to truncate
        t_names <- paste(unname(sapply(names(seq_trunc), crop_names)),
                            '_truncated', sep = '')
        names(seq_trunc) <-
            t_names #new names for sequences to be truncated
        seq_trunc <- subseq(seq_trunc, 1, threshold - 1)
        # shuffle AAStringset to avoid having all the heavy sequences
        # in the last chunk
        seq_set <- sample(c(seq_keep, seq_trunc))

        if (all(width(seq_set) < threshold))
            return(seq_set)
    }
}
