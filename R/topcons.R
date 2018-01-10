#' predict transmembrane domains with TOPCONS
#'
#' This function calls local TOPCONS to predict transmembarne domains in protein
#' sequence. \cr
#' \cr
#' @param input_obj an instance of SignalpResult class, \cr
#'                  input should contain mature_fasta slot;
#' @param TM  allowed number of TM domains in mature peptides,
#' recommended value <= 1; use TM = 0 for strict filtering    
#' @param run_mode piper/starter
#' @param topcons_mode WSDL API / stand alone                
#' @param paths if topcons is not acessible globally, a file
#' conatining a full path to it's executable should be provided; for details
#' please check SecretSanta vignette.
#' @export
#' @return TopconsResult object
#' @examples 

topcons <- function(input_obj,
                    run_mode = c('piper', 'starter'),
                    paths = NULL,
                    topcons_mode = c('API', 'stand-alone'),
                    TM){
    
    # ----- Check the input parameters:
    
    # check that input object belongs to a valid class
    if (is(input_obj, "CBSResult")) {
    } else {
        stop('input_object does not belong to CBSResult class')
    }
    
    # check run_mode value
    
    if (missing(run_mode)) {
        stop('missing argument: run_mode')
    }
    
    run_mode <- match.arg(run_mode)
    
    # check topcons mode
    
    if (missing(topcons_mode)) {
        stop('missing argument: topcons_mode')
    }
    
    topcons_mode <- match.arg(topcons_mode)
    
    # check TM threshold value
    if (!(is.numeric(TM))) stop('TM argument should be numeric')
    if (TM >= 2) warning('Recommended TM threshold values for mature peptides is 0')
    
    # All checked, produce an encouragig message
    message(paste("running topcons in the", topcons_mode, "mode"))
    
    
    #### simple function to run TOPCONS WSDL API script: submit requests to the
    #### web-server
    
    
    #### simple function to run TOPCONS in the stand-alone mode
    
    simple_topcons <- function(fasta) {
        
    }
    
    
    
    
}



