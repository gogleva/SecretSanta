#' run_WolfPsort function
#'
#' This function runs WoLF PSORT to predict protein cellular sub-localisation
#' and returns the most probbale one. Including this step in secretome 
#' prediction pipelines provides additional supportig evidence that a protein
#' might be secreted and deposited outside the cell.\cr
#' \cr
#' Recommended to run on the late stages of secretome prediction pipeline.\cr
#' \cr
#' Also see targetp function - for similar functionality.
#' 
#' @param input_obj Object of CSBResult class
#' @param organism  set relevant taxonomic group,
#'                  options include: \strong{plant},
#'                  \strong{animal}, \strong{fungi};
#' @param paths   if wolfpsort is not acessible globally, a file
#' conatining a full path to it's executable should be provided; for details
#' please check SecretSanta vignette. 
#' @return object of WolfResult class  
#' @export
#' @examples
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata",
#'                                   "sample_prot_100.fasta",
#'                                   package = "SecretSanta"))
#' 
#' # assign this object to the input_fasta slot of CBSResult object
#' inp <- CBSResult(in_fasta = aa[1:10])
#' 
#' # run signalp2 on the initial file:
#' step1_sp2 <- signalp(inp,
#'                      version = 2,
#'                      organism ='euk',
#'                      run_mode = "starter")
#' # run wolfpsort on the signalp output:
#' w <- wolfpsort(step1_sp2, 'fungi')

wolfpsort <- function(input_obj,
                      organism = c('plant', 'animal', 'fungi'),
                      paths = NULL) {
  
  # check that inputs are valid
  
  # check organism argument
  if (missing(organism)) {stop('missing argument: organism')}
  organism <- match.arg(organism)
  
  # check input_obj
  if (is(input_obj, "CBSResult")) {} else {
    stop('input_object does not belong to CBSResult superclass')}
  if (length(getOutfasta(input_obj)) == 0) {
    stop('the input object contains empty out_fasta slot')}

  message("running WoLF PSORT locally...")
  
  fasta <- getOutfasta(input_obj)
  message(paste("Number of submitted sequences...", length(fasta)))

  out_tmp <- tempfile()
  Biostrings::writeXStringSet(fasta, out_tmp)
  
  # get and check paths to wolfpsort
  if (is.null(paths)) {
    full_pa <- 'wolfpsort'
  } else {
    mp <- suppressMessages(manage_paths(in_path = FALSE,
                                        test_mode = 'wolfpsort',
                                        path_file = paths))
    full_pa <- mp$path_tibble$path
  } 
  
  wolf <- system(paste(full_pa,  organism, '<', out_tmp), intern = TRUE)
  
  #parse wolf output
  clean_strings <- function(x, field){
    unlist((stringr::str_split(x, " ")))[c(field)]}
  
  gene_id <- sapply(X = wolf,
                    field = 1,
                    clean_strings,
                    USE.NAMES = FALSE)
  localization <- sapply(X = wolf,
                         field = 2,
                         clean_strings,
                         USE.NAMES = FALSE)
  
  #assemble result tibble with gene id and most probable
  #subsellular localisation
  wolf_tbl <- tibble::as_tibble(data.frame(gene_id, localization)) %>% 
                                  filter_(~gene_id != '#') %>%
                                  filter_(~localization == 'extr')
  
  message(paste('Candidate sequences with extracellular localisation...',
                nrow(wolf_tbl)))

  #assemble wolf result object:
  out_obj <- WolfResult(in_fasta = fasta,
                        out_fasta = fasta[wolf_tbl$gene_id], 
                        wolf_tibble = wolf_tbl)
  
  if (validObject(out_obj)) {return(out_obj)}
}
