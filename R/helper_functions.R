#' crop_names function
#'
#' helper function crop long names for AAStringSet object, return
#' character vector
#' @return character vector
#' @keywords internal
## @export

crop_names <- function(x) { return(unlist(stringr::str_split(x, " "))[1]) }

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
    drop_n <- length(seq_set[Biostrings::width(seq_set) >= threshold])

    if (drop_n == 0)
        return(seq_set)

    if (!truncate) {
        seq_set <- seq_set[Biostrings::width(seq_set) < threshold]
        warning(paste(drop_n, ' long sequenses have been thrown away ...'), call. = FALSE)
        return(seq_set)

    } else if (truncate) {
        message(paste(drop_n, ' sequences need to be truncated ...'))
        seq_keep <-
            seq_set[Biostrings::width(seq_set) < threshold] # not so long sequences
        seq_trunc <-
            seq_set[Biostrings::width(seq_set) >= threshold] # sequences to truncate
        t_names <- paste(unname(sapply(names(seq_trunc), crop_names)),
                            '_truncated', sep = '')
        names(seq_trunc) <-
            t_names #new names for sequences to be truncated
        seq_trunc <- subseq(seq_trunc, 1, threshold - 1)
        # shuffle AAStringset to avoid having all the heavy sequences
        # in the last chunk
        seq_set <- sample(c(seq_keep, seq_trunc))

        if (all(Biostrings::width(seq_set) < threshold))
            return(seq_set)
    }
}

#' estimate_lim
#'
#' helper function to estimate approximate length threshold if
#' sequence chunk size exceedes 200000 a.a len limit and truncate long sequences
#' or throw them away, otherwise signalp will break (signalp2 and signalp3 will)
#' @return truncated AAStringSet
#' @param fasta_chunk AAStringSet
#' @param truncate truncate or throw away
#' @keywords internal
## @export

estimate_lim <- function(fasta_chunk, truncate) {

    len_lim <- (200000 / length(fasta_chunk) + 50)
    if (sum(Biostrings::width(fasta_chunk)) >= 200000) {
        message(paste(
                'fasta size exceedes maximal total residue limit, seqs > ',
                round(len_lim),
                ' residues will be truncated ...'
            )
        )
        fasta_trunc <-
            truncate_seq(truncate = truncate, fasta_chunk, len_lim)
        return(fasta_trunc)
    } else {
        return(fasta_chunk)
    }
}
