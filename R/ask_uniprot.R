#' fetch known subcellular location from UniprotKB based on uniprot ids
#'
#' @param uniprotID_list  a list of uniprot IDs
#' @export 
#' @return a list of locations
#' @examples 
#' id_list <- c('P39864', 'D0N4E2', 'Q5BUB4', 'D0N381', 'B1NNT7', 'D0NP26')
#' res <- ask_uniprot(id_list)
#' #try submitting a list containing some non-existing UniprotIDs:
#' bad_list <- c('P39864', 'D0N4E2', 'Q5BUB4', 'D0N381', 'B1NNT7', 'D0NP2688', 'D0N4E2222')
#' res2 <- ask_uniprot(bad_list)

ask_uniprot <- function(uniprotID_list) {
    
    # convert long ids to short uniprot ids:
    mnchar <- mean(sapply(uniprotID_list, nchar, simplify = TRUE)) > 10
    pipech <- sapply(uniprotID_list, stringr::str_count, pattern = "\\|",
                     USE.NAMES = FALSE) 
        
    if (mnchar & all(pipech > 1)) { 
        # extract short uniprot ids
        message('Detected long ids, converting to short UniprotIDs ...')
        uniprotID_list <- sapply(uniprotID_list, function(x) str_split(x, "\\|")[[1]][2], simplify = TRUE, USE.NAMES = FALSE)
    }

    simple_fetch <- function(uniprotID) {
        base <- 'http://www.uniprot.org/uniprot/'
        fetch_url <- paste(base, uniprotID, '.txt', sep = '')
        resp <- httr::GET(fetch_url)
        
        # check here for status errors:
        if (httr::http_error(resp)) {
            msg <- httr::http_status(resp)$message
            keys <- loc <- go_cc <- msg
            warning(paste('Location retrival for provided UniprotID', uniprotID,
                          'failed due to:', '\n', msg,
                          '\n', 
                          'please check the offending id.', sep = ' '), 
                    call. = FALSE)
            
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
                keys <- stringr::str_replace(raw_keys, kw_pattern, '') %>%
                    paste(collapse = ' ')
                
            } else {
                keys <- 'not present'
            }
            
            # fetch relevant GO terms (CC)
            go_pattern <- '^DR   GO; '
            
            if (any(stringr::str_detect(cont, go_pattern))) {
                go <- cont[stringr::str_detect(cont, go_pattern)]
                go_cc <- stringr::str_replace(go, go_pattern, '')
                go_cc <- go_cc[grepl('C:', go_cc)] %>%
                    paste(collapse = ' ')
                
            } else {
                go_cc <- 'not present'
            }
        }
        
        result <- tibble(UniprotID = uniprotID,
                         Subcellular.Location = loc,
                         Key.Words = keys,
                         GO.CC = go_cc)
    }
    
    message("Fetching location from UniprotKB ...")
    all_res <- dplyr::bind_rows(lapply(uniprotID_list, simple_fetch))
}
