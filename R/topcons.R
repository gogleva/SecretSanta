#' parse TOPCONS output with prediction of transmenbrane domains
#'
#' This function parses TOPCONS outputs.
#' Will be restricted to work in the piper mode \cr
#' \cr
#' @param input_obj an instance of CBSResult class, produced by previous step;
#' @param TM  allowed number of TM domains in mature peptides,
#' recommended value <= 1; use TM = 0 for strict filtering    
#' @param parse_dir dir with archived topcons output
#' @param topcons_mode WEB / WSDL-API / stand-alone              
#' @export
#' @return TopconsResult object
#' @examples 

topcons <- function(input_obj = NULL,
                    parse_dir,
                    topcons_mode = c('API', 'WEB-server', 'stand-alone'),
                    TM){
    
    # ----- Check the input parameters:
    
    # check run_mode value
    
    if (missing(run_mode)) {
        stop('missing argument: run_mode')
    }
    
    run_mode <- match.arg(run_mode)
    
    # check that path to zipped output is provided:
    
    if (missing(parse_dir) {
        stop('missing argument: topcons_mode')
    }
    
    # check that topcons output file exists:
    
    if (!file.exists(parse_dir)) {
        stop('Please provide valid path to the zipped TOPCONS output')
    }   
    
    # check that input object belongs to a valid class,
    
    # expect CBSResult object    
    if (!(is(input_obj, "CBSResult"))){
        stop('input_object does not belong to CBSResult class')
    } 

    # check topcons mode
    
    if (missing(topcons_mode)) {
        stop('missing argument: topcons_mode')
    }
    
    topcons_mode <- match.arg(topcons_mode)
    
    # check TM threshold value
    if (!(is.numeric(TM))) stop('TM argument should be numeric')
    if (TM >= 2) warning('Recommended TM threshold values for mature peptides is 0')
    
    # check that input_obj contains non-empty out_fasta slot
   if (run_mode == 'piper') {
        if (length(getOutfasta(input_obj)) != 0) {
            fasta <- getOutfasta(input_obj)
        } else {
            stop('out_fasta attribute is empty')
        }
    }
    
    
    # All checked, produce an encouragig message
    message(paste("running topcons in the", topcons_mode, "mode"))

    
    parse_topcons <- function(dir_to_parse) {
        
        
        #----- Run topcons prediction:
        message(paste('Number of submitted sequences...', length(aaSet)))
        
        # convert fasta to temp file:
        out_tmp <- tempfile()
        Biostrings::writeXStringSet(aaSet, out_tmp) #write tmp fasta file
        
        
        if (topcons_mode == 'API') {
            # we need topcons WSDL API script
            
            # get and check paths for topcons WSDL API script
            # submit requests to the web-server, this option is not 
            # very fast
        } else if {topcons_mode == 'stand-alone'}

        
        
        
        
    #### run TOPCONS in the stand-alone mode    
    
        
    }

}



