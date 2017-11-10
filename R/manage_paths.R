#' organise and test external dependencies
#'
#' This function checks that all the external dependencies required for
#' secertome prediction pipelines are available and can sucessfully pass test
#' runs. \cr
#' \cr
#' All the external dependencies can be made accessible via $PATH environment
#' variable; alternatively a text file the with full paths to all the external
#' dependencies should be provided.
#' \cr
#' \cr
#' For the instalation instructions please see SecretSanta vignette.
#' @param in_path if \strong{TRUE} manage_paths will attempt to run
#' external dependencies, assuming that respective paths are attached to the
#' $PATH variable;\cr
#' \cr
#' if \strong{FALSE} you should supply path_file.
#' @param test_tool if \strong{all} - all the external dependencies will be
#' checked;\cr
#' alternatively specify a tool name to be checked.
#' @param path_file    full paths to external dependencies in a
#' 2-column space-separated text file;\cr
#' for multiple versions of signalp please use 'signalpV'notation, where V is a
#' version number;\cr
#' \cr
#' \strong{first column} should contain tool name; \cr
#' \strong{second column} should contain full path to the tool's executable.
#' @return a list of length 3 with the following elements:
#' \itemize{
#' \item \strong{tests}    TRUE if all the exteranl dependencies are working;
#' \item \strong{in_path}    TRUE if the dependecies are accessible via $PATH;
#' FALSE - if paths are provided with the path_file argument;
#' \item \strong{path_tibble} a tibble with verified names and paths for
#' external dependencies; \cr
#' NA if in_path is TRUE;\cr
#' if a specific tool is tested, a singualr path will be returned.
#' }
#' @export
#' @examples
#' # Example1:
#' # here we assume that paths to all the
#' # external dependencies are attached to
#' # the $PATH variable:
#'
#' manage_paths(in_path = TRUE, test_tool = 'all')
#'
#' # to test just a single tool:
#' manage_paths(in_path = TRUE, test_tool = 'signalp2')
#'
#' # Example2:
#' # alternatively, we are in a situation
#' # when changing $PATH is not possible,
#' # so we supply a file with listed full
#' # paths to the external dependencies:
#'
#' manage_paths(in_path = FALSE,
#' test_tool = 'all',
#' path_file = system.file("extdata", "sample_paths", package = "SecretSanta"))

manage_paths <- function(in_path = c(TRUE, FALSE),
                            test_tool = c('all', 'signalp2', 'signalp3',
                                            'signalp4','targetp', 'tmhmm',
                                            'wolfpsort'),
                            path_file = NULL) {
    
    test_tool <- match.arg(test_tool)
    
    # here we assume that all the tools are accessible via $PATH:
    if (in_path == TRUE) {
        message('checking dependencies accessible via $PATH')
        if (is.null(path_file)) {
        } else {
            message('path file provided, but not required')
        }
        pp <- NA
    }
    
    # alternatively, we will be dealing with paths supplied in a separate
    # path_file, first check that all the supplied paths are valid
    if (in_path == FALSE) {
        if (is.null(path_file)) {
            stop('Please supply path_file')
        } else {
            # read path file in a tibble
            pp <- tibble::as.tibble(read.table(path_file, sep = ' ',
                                               header = FALSE,
                                               stringsAsFactors = FALSE))

            # check that there are only 2 columns, i.e no extra spaces in
            # tool names
            if (!(ncol(pp) == 2)) {
                stop('Please check that there are no spaces in the tool names.')
            }
            
            names(pp) <- c("tool", "path")
            # check that all supplied paths exist
            pp$status <- file.exists(pp$path)
            
            if (all(pp$status)) {
                message('All paths are valid')
                # convert all the tool names to lower case to avoid confusion
                pp$tool <- tolower(pp$tool)
                
            } else {
                message('Supplied file path does not exist.')
                message(sapply(pp[pp$status == TRUE,]$path, paste, '\n'))
                stop('Please check that supplied paths are correct')
            }
        }
        pp
    }
    
    # check that all the dependencies are executable in principle,
    # i.e we are able to process a small sample fasta file
    
    # micro fasta file to test with all tools:
    test_fasta <- system.file("extdata", "small_prot.fasta", 
                                package = "SecretSanta")
    
    # helper function to extract tool paths when provided in path_file

    get_paths <- function(tool_name) { pp[pp$tool == tool_name,]$path }
    
    # helper function to generate status message: tool failed or not

    status_message <- function(tool_name, status) {
        if (status == 'OK') {message(paste(tool_name, 'run completed'))}
        if (status == 'FAIL') {message(paste(tool_name, 'test_run_failed'))}
    }
        
    # helper function to make a tool call, wrap in get_paths if necessary
    make_call <- function(tool) {if (in_path == TRUE) tool else get_paths(tool)}
        
    # now we will wrap calls and evaluation of the outputs in small functions:
    
    ## need to re-write this as a list of functions?/function factory?
    ## test_fasta will be in parent function environment, so we may skip
    ## this argument
    
    test_my_tool <- function(tool_name, call_param, grep_param, grep_line){
        tool_call <- system(paste(make_call(tool_name),
                                  call_param,
                                  test_fasta), intern = TRUE)
        if (grepl(grep_param, tool_call[grep_line]))  {
            status_message(tool_name, 'OK')
            tool_status <- TRUE
        } else {
            status_message(tool_name, 'FAIL')
            tool_status <- FALSE
        }
        return(tool_status)
    }

    # gather parameters in lists
    tool_names <- c('signalp2', 'signalp3', 'signalp4', 'targetp', 'tmhmm',
                    'wolfpsort')
    call_params <- c(rep('-t euk', 3), '-P', '--short', 'fungi <')
    grep_params <- c('SignalP predictions',
                     'SignalP 3.0 predictions',
                     'SignalP-4',
                     'targetp v1.1 prediction results',
                     'ALI_PLTG_1\tlen=94\tExpAA=22.44',
                     '# k used for kNN is: 27')
    grep_lines <- c(rep(1, 3), 2, rep(1,2))

    # now we will run required tests according to the actual value
    # of the mode argument
    
    if (test_tool == 'all') {
        all_tests <- Map(test_my_tool, tool_names, call_params,
                       grep_params, grep_lines)
    } else {
        t_ind <- tool_names == test_tool
        all_tests <- test_my_tool(test_tool, call_params[t_ind],
                             grep_params[t_ind],
                             grep_lines[t_ind])
    } 
    
    # construct the final output
    
    # # output paths only for the tools tested
    # if (in_path == FALSE & (test_tool != 'all')) {
    #     pp <- pp[pp$tool == test_tool,] 
    # }
    
#    result <- list(tests = all_tests, in_path = in_path, path_tibble = pp)
    
    return(all_tests)
}
