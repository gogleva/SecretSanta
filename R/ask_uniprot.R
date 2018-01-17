#' fetch known subcellular location from UniprotKB based on uniprot ids
#'
#' @param uniprotID_list  a list of uniprot IDs
#' @export 
#' @return a list of locations
#' @examples 
#' id_list <- c('P39864', 'D0N4E2', 'Q5BUB4', 'D0N381', 'B1NNT7', 'D0NP26')
#' res2 <- ask_uniprot(id_list)
#' #try submitting a list with non-existing UniprotIDs:
#' bad_list <- c('P39864', 'D0N4E2', 'Q5BUB4', 'D0N381', 'B1NNT7', 'D0NP2688', 'D0N4E2222')
#' res2 <- ask_uniprot(bad_list)

ask_uniprot <- function(uniprotID_list) {
    
    simple_fetch <- function(uniprotID) {
        base <- 'http://www.uniprot.org/uniprot/'
        fetch_url <- paste(base, uniprotID, '.txt', sep = '')
        resp <- httr::GET(fetch_url)
        
        # check here for status errors:
        if (httr::http_error(resp)) {
            msg <- httr::http_status(resp)$message
            keys <- loc <- msg
            warning(paste('location retrival for provided UniprotID', uniprotID,
                          'failed due to:', '\n', msg,
                          '\n', 
                          'please check the offending id', sep = ' '), call. = FALSE)
            
        } else {
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
        }
        
        result <- tibble(UniprotID = uniprotID,
                         Subcellular.Location = loc,
                         Key.Words = keys)
    }
    
    all_res <- dplyr::bind_rows(lapply(uniprotID_list, simple_fetch))
}
