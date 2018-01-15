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

topcons <- function(input_obj,
                    parse_dir,
                    topcons_mode = c('API', 'WEB-server', 'stand-alone'),
                    TM,
                    SP = FALSE){
    
    # ----- Check the input parameters:
    
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
        'Recommended TM threshold values for secreted peptides is 0')
    
    # check that input_obj contains non-empty out_fasta slot
    if (length(getOutfasta(input_obj)) != 0) {
            fasta <- getOutfasta(input_obj)
        } else {
            stop('out_fasta attribute is empty')
        }

    # All checked, produce an encouragig message
    message(paste("running topcons parser for", topcons_mode, "format"))
    
    # for tests:
    dir_to_parse <- "/home/anna/anna/Labjournal/SecretSanta_external/rst_SVw4hG.zip"
    fasta <- getOutfasta(sp)

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
        
        # remove unpacked dir
        unlink(paste(dirname(dir_to_parse), '/', rst_id, sep = ''),
               recursive = TRUE)
        
        names(topcons_tibble) <- c('seq', 'length', 'TM',
                                   'SP', 'source', 'run_time',
                                   'gene_id')
        
        # crop names, to remove additional annotations:
        topcons_tibble$gene_id <- unname(sapply(topcons_tibble$gene_id, crop_names))
        
        # keep gene_ids only present in the input object
        topcons_tibble <- topcons_tibble[topcons_tibble$gene_id %in% names(fasta),]
        
        # filter based on TM threshold
        topcons_tibble <- (topcons_tibble %>% dplyr::filter_( ~ TM <= TM))
        
        # filter based on SP threshold
        if (SP == TRUE) {
        topcons_tibble <- topcons_tibble %>% dplyr::filter_( ~ SP == 'True')
        }
        
        # assemble TOPCONS object
        
        out_obj <- TopconsResult(top_tibble = topcons_tibble)
        out_obj <- setInfasta(out_obj, in_fasta = fasta)
        out_obj <- setOutfasta(out_obj, out_fasta = fasta[topcons_tibble$gene_id])

        if (validObject(out_obj)) { return(out_obj)}

    }
    
    parse_topcons(parse_dir)

}


# tests: 
#dir_to_parse = "/home/anna/anna/Labjournal/SecretSanta_external/TOPCONS2_API/rst_ArIQg1.zip"







