#' fetch known subcellular location from UniprotKB based on uniprot ids
#'
#' @param UniprotID_list  a list of UniprotIDs
#' @export 
#' @return a list of locations
#' @examples 
#' id_list <- c('P39864', 'D0N4E2', 'Q5BUB4', 'D0N381', 'B1NNT7', 'D0NP26')
#' known_locations <- ask_uniprot(id_list)

ask_uniprot <- function(uniprotID_list) {
    
    simple_fetch <- function(uniprotID) {
        base <- 'http://www.uniprot.org/uniprot/'
        fetch_url <- paste(base, uniprotID, '.txt', sep = '')
        resp <- httr::GET(fetch_url)
        cont <- stringr::str_split(httr::content(resp, "text"),
                                   '\n',
                                   simplify = TRUE)
        
        # fetch subcellular location field
        loc_pattern <- 'SUBCELLULAR LOCATION:'
        if (any(stringr::str_detect(cont, loc_pattern))) {
            raw_loc <- cont[stringr::str_detect(cont, loc_pattern)]
            loc <- stringr::str_replace(raw_loc,
                                        'CC   -!- SUBCELLULAR LOCATION: ', '')
        } else {
            loc <- 'not present'
        }
        
        # fetch and parse keywords
        kw_pattern <- '^KW   '
        if (any(stringr::str_detect(cont, kw_pattern))) {
            raw_keys <- cont[stringr::str_detect(cont, kw_pattern)]
            keys <- stringr::str_replace(raw_keys, 'KW   ', '') %>%
                paste(collapse = ' ')
            
        } else {
            keys <- 'not present'
        }
        
        result <- list(loc, keys)
    }
    
    all_res <- sapply(uniprotID_list, simple_fetch, simplify = TRUE)
}

### experiments
