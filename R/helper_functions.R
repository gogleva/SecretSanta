#' crop_names function
#'
#' helper function crop long names for AAStringSet object, return
#' character vector
#' @return character vector 
#' @export 
#' @keywords internal

crop_names <- function(x){unlist(stringr::str_split(x, " "))[1]}


