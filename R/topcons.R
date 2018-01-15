#' parse TOPCONS output with prediction of transmenbrane domains
#'
#' This function parses TOPCONS outputs.
#' Will be restricted to work in the piper mode \cr
#' \cr
#' @param input_obj an instance of CBSResult class, produced by previous step;
#' @param TM  allowed number of TM domains in mature peptides,
#' recommended value <= 1; use TM = 0 for strict filtering
#' @param SP do we need to filter for presence of SP? probably yes? By default - TRUE
#' @param parse_dir dir with archived topcons output
#' @param topcons_mode WEB / WSDL-API / stand-alone              
#' @export
#' @return TopconsResult object
#' @examples 

topcons <- function(input_obj = NULL,
                    parse_dir,
                    topcons_mode = c('API', 'WEB-server', 'stand-alone'),
                    TM,
                    SP = TRUE){
    
    # ----- Check the input parameters:
    
    # check run_mode value
    
    if (missing(run_mode)) { stop('missing argument: run_mode')}
    run_mode <- match.arg(run_mode)
    
    # check that path to zipped output is provided:
    
    if (missing(parse_dir)) {stop('missing argument: topcons_mode')}
    
    # check that topcons output file exists:
    
    if (!file.exists(parse_dir)) {
        stop('Please provide valid path to the zipped TOPCONS output')
    }   
    
    # check that input object belongs to a valid class,
    if (!(is(input_obj, "CBSResult"))){
        stop('input_object does not belong to CBSResult class')
    } 

    # check topcons mode
    if (missing(topcons_mode)) {stop('missing argument: topcons_mode')}
    topcons_mode <- match.arg(topcons_mode)
    
    # check TM threshold value
    if (!(is.numeric(TM))) stop('TM argument should be numeric')
    if (TM >= 2) warning(
        'Recommended TM threshold values for mature peptides is 0')
    
    # check that input_obj contains non-empty out_fasta slot
    if (length(getOutfasta(input_obj)) != 0) {
            fasta <- getOutfasta(input_obj)
        } else {
            stop('out_fasta attribute is empty')
        }

    # All checked, produce an encouragig message
    message(paste("running topcons in the", topcons_mode, "mode"))

    parse_topcons <- function(dir_to_parse) {
        
        # first, unzip the archive
        rst_id <- strsplit(basename(dir_to_parse), split = '.zip')[[1]]
        unzp <- unzip(zipfile = dir_to_parse, 
                      exdir = paste(dirname(dir_to_parse),
                                    '/',
                                    rst_id,
                                    sep = ''))
        # filter for file name with summary stats:
        unzp <- unzp[grepl('finished_seqs.txt', unzp)]
        topcons_tibble <- tibble::as.tibble(read.table(unzp,
                                                       sep = '\t',
                                                       header = FALSE,
                                                       quote = ""))
        names(topcons_tibble) <- c('seq', 'length', 'TM',
                                   'SP', 'source', 'run_time',
                                   'seqID')
        
        # crop names, to remove additional annotations:
        topcons_tibble$seqID <- sapply(crop_names, topcons_tibble$seqID)
        
        # filter based on TM threshold
        
        
        # filter based on SP threshold
        
        
        # assemble TOPCONS object
        
    
    }

}


# tests: 
#dir_to_parse = "/home/anna/anna/Labjournal/SecretSanta_external/TOPCONS2_API/rst_ArIQg1.zip"

#unzip(zipfile = "/home/anna/anna/Labjournal/SecretSanta_external/TOPCONS2_API/rst_ArIQg1.zip", exdir= di"/home/anna/anna/Labjournal/SecretSanta_external/TOPCONS2_API/rst_unzip")








