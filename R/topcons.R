# helper function to extract number of TM domains or check for presence of SP 
# based on stand-alone TOPCONS predictions
# sample input string: ('iiioooMMMMMMMMMMMMMMMMooooMMMMMMMiiioooo')

stretch_parse <- function(str, core_pattern){
    runs <- paste(rle(strsplit(str, "")[[1]])$values, collapse="")
    str_count(runs, core_pattern)
}

#' parse TOPCONS output 
#' 
#' This function parses results of TOPCONS method, used to predict transmembrane domains. The parser is restricted to 'piper'-like behaviour.\cr
#' @param input_obj an instance of CBSResult class, produced by previous step;
#' @param TM  allowed number of TM domains in mature peptides,
#' recommended value <= 1; use TM = 0 for strict filtering
#' @param SP filter according to TOPCONS prediction of signal peptides. TRUE - keep only proteins containg signal peptides according to TOPCONS; FALSE - disable filtering for TOPCONS-based predictions of signal peptides. Deafault = FALSE.
#' @param parse_dir dir with archived topcons output (for WSDL-API or WEB outputs)
#' or dir with output produced by stand-alone version.
#' @param topcons_mode
#' \itemize{
#' \item    WEB - output from TOPCOS2 web-server \url{http://topcons.net/}
#' \item    WSDL-API - output returned by WSDL API script \url{http://topcons.net/pred/help-wsdl-api/}
#' \item    stand-alone - results of TOPCONS stand-alone version \url{https://github.com/ElofssonLab/TOPCONS2}
#' }
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
    
    if ((topcons_mode %in% c('API', 'WEB-server')) & (!file.exists(parse_dir))) {
        stop('Please provide valid path to the zipped TOPCONS output')
    }
    
    if (topcons_mode == 'stand-alone') {
        if (!(file.exists(paste(parse_dir, 'time.txt', sep = '')))) {
            stop('could not find "time.txt" file, please check the path')
        }
        if (!(any(grepl('seq_', list.files(parse_dir))))) {
            stop('could not find "seq_*" subdirectories, please check the path')
        }
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
    
    if (!(is.logical(SP))) {stop('SP argument must be logical')}

    # All checked, produce an encouragig/status messages
    
    message(paste("running topcons parser for", topcons_mode, "format"))
    
    if (SP == FALSE) {
        message('verify SP predictions ... YES')
    } else {
        message("verify SP predictions ... NO")
    }
       
    message(paste('TM domains allowed ... '), TM)
    
    #essential fucntion to run simple topcons parser:
    
    parse_topcons <- function(dir_to_parse, topcons_mode) {
        
        # parser bifurcation based on input format:
        
        if (topcons_mode %in% c("API", "WEB-server")) {
        
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
        topcons_tibble$gene_id <- sapply(topcons_tibble$gene_id, crop_names,
                                         USE.NAMES = FALSE)
        }
        
        # stand-alone mode does not produce summary table
        # we need to collate/extract relevant fields from isolated files and 
        # dirs to match WEB/API format
        
        if (topcons_mode == 'stand-alone') {

            # extract seq_ids and run times from time.txt
            t_file <- read.table(paste(dir_to_parse, '/time.txt', sep = ''),
                                 sep = ';')
            # crop names, just in case
            gene_id <- sapply(t_file$V1, crop_names)
            run_times <- t_file$V2
            
            # extract lengths from seq_*/dg.txt files -> combine
            
            dirs <- list.files(path = dir_to_parse, pattern = "seq_")
            
            # make sure that ordering of seq_folders is correct, we don't need the
            # default ordering
            vals <- as.numeric(gsub("seq_", "", dirs))
            dirs <- dirs[order(vals)]
            
            # exact paths to dg.txt files
            dgs  <- sapply(dirs, function(x) paste(dir_to_parse,
                                                   x,
                                                   'dg.txt',
                                                   sep = '/'))
            
            # helper function to extract length field w/o reading the whole file
            get_length <- function(dg_path) {
                system(paste('cat', dg_path, "| grep SeqLength: | awk '{print $2}'"),
                       intern = TRUE) %>% as.numeric()
            }
            
            lengths <- sapply(dgs, get_length)
                    

            # extract number of TM and SP domains
            
            tops <- sapply(dirs, function(x) paste(dir_to_parse,
                                                        x,
                                                        'Topcons/topcons.top',
                                                        sep = '/'))
                
            # helper function to extract raw preds in correct order:
            get_preds <- function(top_path) {
                system(paste("cat", top_path, sep = ' '), intern = TRUE)
            }
            
            raw_preds <- sapply(tops, get_preds)
            
            # extract number of TM domains:
            TM_num <- sapply(raw_preds, function(x) stretch_parse(x, 'M'),
                         USE.NAMES = FALSE)
            
            SP_if <- sapply(raw_preds, function(x) stretch_parse(x, 'S'),
                         USE.NAMES = FALSE) %>% as.logical()
            
            # assemble the tibble now
            topcons_tibble <- tibble::tibble('seq' = dirs,
                                             'length' = lengths,
                                             'TM' = TM_num,
                                             'SP' = SP_if,
                                             'source' = 'stand-alone',
                                             'run_time' = run_times,
                                             'gene_id' = gene_id
                                             )
            
        }
    
        # ----filter the tibble and assemble result object
    
        # keep gene_ids only present in the input object
        topcons_tibble <- topcons_tibble[topcons_tibble$gene_id %in% names(fasta),]
    
        # filter based on TM threshold
        topcons_tibble <- (topcons_tibble %>% dplyr::filter_( ~ TM <= TM))
    
        # filter based on SP threshold
        if (SP == TRUE) {
            topcons_tibble <- topcons_tibble %>% dplyr::filter_( ~ SP == 'True')
        }
    
        message(paste('Number of result proteins ...', nrow(topcons_tibble)))
        # assemble TOPCONS object
    
        out_obj <- TopconsResult(top_tibble = topcons_tibble)
        out_obj <- setInfasta(out_obj, in_fasta = fasta)
        out_obj <- setOutfasta(out_obj, out_fasta = fasta[topcons_tibble$gene_id])
    
        if (validObject(out_obj)) { return(out_obj)}

    }
    
    parse_topcons(parse_dir)

}

#dir_to_parse <- "/home/anna/anna/Labjournal/SecretSanta_external/TOPCONS2_stand-alone/rst_milti/multiple_seqs"
#inp <- readAAStringSet("")
    
    
    
