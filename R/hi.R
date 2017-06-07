#' Hi function
#'
#' This is a test fucntion
#' @param name enter your user name here
#' @keywords hi
#' @export
#' @examples
#' hi()

hi <- function(user_name) {
  print(paste("hi", user_name, "let's get started!", sep = " "))
  }

#message for the user generate at the moment of package build
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to SecretSanta!")
}

