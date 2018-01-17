#' fetch known subcellular location from UniprotKB based on uniprot ids
#'
#' @param UniprotID_list  a list of UniprotIDs
#' @export 
#' @return a list of locations
#' @examples 
#' id_list <- c('P39864', 'D0N4E2', 'Q5BUB4', 'D0N381', 'B1NNT7', 'D0NP26')
#' known_locations <- fetch_uniprot_loc(id_list)

fetch_uniprot_loc <- function(uniprotID_list) {
    
    simple_fetch <- function(uniprotID) {
        base <- 'http://www.uniprot.org/uniprot/'
        fetch_url <- paste(base, uniprotID, '.txt', sep = '')
        resp <- httr::GET(fetch_url)
        cont <- stringr::str_split(content(resp, "text"), '\n', simplify = TRUE)
        match_pattern <- 'SUBCELLULAR LOCATION:'
        if (any(stringr::str_detect(cont, match_pattern))) {
            raw_loc <- cont[stringr::str_detect(cont, match_pattern)]
        } else {
            raw_loc <- 'not present'
        }
        raw_loc
        }
    
    all_res <- sapply(uniprotID_list, simple_fetch, simplify = TRUE)
}
