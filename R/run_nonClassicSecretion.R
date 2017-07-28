#' Prediction of non-classical protein secretion
#'
#' This function runs SecretomeP 2.0 to produce ab initio predictions of non-classical i.e. not signal peptide triggered protein secretion. 
#' 
#' @param input_obj input object of CBSResult superclass
#' @export run_mode 'starter' or 'piper' 
#' @examples 
#' 

nonClassicSecretion <- function(input_obj, run_mode) {
  # check the input
  if (is(input_obj, "CBSResult")) {} else stop('input_object does not belong to CBSResult superclass')
}
