#' manage_paths function
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
#' @param test_mode if \strong{all} - all the external dependencies will be
#'                checked;\cr
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
#' \item \strong{in_path}    TRUE if the dependecies are acessible via $PATH;
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
#' manage_paths(in_path = TRUE, test_mode = 'all')
#'
#' # to test just a single tool:
#' manage_paths(in_path = TRUE, test_mode = 'signalp2')
#'
#' # Example2:
#' # alternatively, we are in a situation
#' # when changing $PATH is not possible,
#' # so we supply a file with listed full
#' # paths to the external dependencies:
#'
#' manage_paths(in_path = FALSE,
#' test_mode = 'all',
#' path_file = system.file("extdata",
#'              "sample_paths",
#'              package = "SecretSanta"))

manage_paths <- function(in_path = c(TRUE, FALSE),
                        test_mode = c('all','signalp2','signalp3',
                                      'signalp4', 'targetp', 'tmhmm',
                                      'wolfpsort'),
                        path_file = NULL) {
    test_mode <- match.arg(test_mode)
    
    # here we assume that all the tools are acessible via $PATH:
    if (in_path == TRUE) {
        message('checking dependencies acessible via $PATH')
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
            pp <- suppressMessages(readr::read_delim(path_file,                                                                                     delim = ' ',                                                                                      col_names = FALSE))
            # check that there are only 2 columns, i.e no extra spaces in
            # tool names
            if (!(ncol(pp) == 2)) {
                stop('Please ensure that there are no spaces in the tool names.')
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
                message(sapply(pp %>%
                        dplyr::filter_( ~status == FALSE) %>%
                        dplyr::select_( ~path),
                    paste, '\n'
                ))
                stop('Please check that supplied paths are correct')
            }
        }
        pp
    }
    
    # check that all the dependencies are executable in principle,
    # i.e we are able to process a small sample fasta file
    
    # micro fasta file to test with all tools:
    test_fasta <- system.file("extdata", "small_prot.fasta",                                                        package = "SecretSanta")
    
    # helper function to extract tool paths when provided in path_file
    get_paths <- function(tool_name) {
        pp %>%
        filter_(~ tool == tool_name) %>%
        select_(~ path)
    }
    
    # helper function to generate success message
    ok_message <- function(tool_name) {
        message(paste(tool_name,
                                    'run completed'))
    }
    
    # helper function to generate failure message
    failure_message <- function(tool_name) {
        message(paste(tool_name, 'test run failed'))
    }
    
    # helper function to make a tool call, wrap in get_paths if necessary
    make_call <- function(tool) {
        if (in_path == TRUE) {
            tool <- tool
        } else if (in_path == FALSE) {
            tool <- get_paths(tool)
        }
    }
    
    # now we will wrap calls and evaluation of the outputs in small functions:
    
    #signalp2:
    test_signalp2 <- function() {
        sp2_call <-
            system(paste(make_call('signalp2'), '-t euk', test_fasta),
                   intern = TRUE)
        if (grepl('SignalP predictions', sp2_call[1])) {
            ok_message('signalp2')
            sp2_test <- TRUE
        } else {
            failure_message('signalp2')
            sp2_test <- FALSE
        }
    }
    
    #signalp3:
    test_signalp3 <- function() {
        sp3_call <-
            system(paste(make_call('signalp3'), '-t euk', test_fasta),
                         intern = TRUE)
        if (grepl('SignalP 3.0 predictions', sp3_call[1])) {
            ok_message('signalp3')
            sp3_test <- TRUE
        } else {
            failure_message('signalp3')
            sp3_test <- FALSE
        }
    }
    
    #signalp4:
    test_signalp4 <- function() {
        sp4_call <-
            system(paste(make_call('signalp4'), '-t euk', test_fasta),
                         intern = TRUE)
        if (grepl('SignalP-4', sp4_call[1])) {
            ok_message('signalp4')
            sp4_test <- TRUE
        } else {
            failure_message('signalp4')
            sp4_test <- FALSE
        }
    }
    
    #targetp:
    test_targetp <- function() {
        tp_call <- system(paste(make_call('targetp'), '-P', test_fasta),
                                            intern = TRUE)
        if (grepl('targetp v1.1 prediction results', tp_call[2])) {
            ok_message('targetp')
            tp_test <- TRUE
        } else {
            failure_message('targetp')
            tp_test <- FALSE
        }
    }
    
    # tmhmm:
    test_tmhmm <- function() {
        tm_call <- system(paste(make_call('tmhmm'), '--short', test_fasta),
                                            intern = TRUE)
        if (grepl('ALI_PLTG_1\tlen=94\tExpAA=22.44', tm_call[1])) {
            ok_message('tmhmm')
            tm_test <- TRUE
        } else {
            failure_message('tmhmm')
            tm_test <- FALSE
        }
    }
    
    # wolfpsort:
    test_wolfpsort <- function() {
        wolf_call <- system(paste(make_call('wolfpsort'), 'fungi', '<',
                                                            test_fasta),
                                                intern = TRUE)
        
        if (grepl('# k used for kNN is: 27', wolf_call[1])) {
            ok_message('wolfpsort')
            wolf_test <- TRUE
        } else {
            failure_message('wolfpsort')
            wolf_test <- FALSE
        }
    }
    
    # now we will run required tests according to the actual value
    # of the mode argument
    
    if (test_mode == 'all') {
        all_tests <- all(
            test_signalp2(),
            test_signalp3(),
            test_signalp4(),
            test_tmhmm(),
            test_targetp(),
            test_wolfpsort()
        )
    } else if (test_mode == 'signalp2') {
        all_tests <- test_signalp2()
    } else if (test_mode == 'signalp3') {
        all_tests <- test_signalp3()
    } else if (test_mode == 'signalp4') {
        all_tests <- test_signalp4()
    } else if (test_mode == 'tmhmm') {
        all_tests <- test_tmhmm()
    } else if (test_mode == 'targetp') {
        all_tests <- test_targetp()
    } else if (test_mode == 'wolfpsort') {
        all_tests <- test_wolfpsort()
    }
    
    # construct the final output
    
    # output paths only for the tools tested
    if (in_path == FALSE & (test_mode != 'all')) {
        pp <- pp %>% filter_( ~ tool == test_mode)
    }
    
    result <- list(tests = all_tests,
                   in_path = in_path,
                   path_tibble = pp)
    
    return(result)
}
