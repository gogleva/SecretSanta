#' run_WolfPsort function
#'
#' This function runs WoLF PSORT to predict protein cellular sub-localisation, returns the most probbale one. Provides additional supportig evidence that a protein might be secreted and deposited outside the cell. Recommended to run on the late stages of secretome prediction pipeline.
#' @param input_obj Object of CSBResult class
#' @param organism  set relevant taxonomic group,
#'                  Options include: plant, animal, fungi;
#' @param paths   tibble with paths to external dependencies, generated with \code{\link{manage_paths}} function
#' @return object of WolfResult class  
#' @export
#' @examples
#' my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
#' 
#' # initialise SignalpResult object
#' inp <- SignalpResult()
#' 
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), use.names = TRUE)
#' 
#' # assign this object to the input_fasta slot of SignalpResult object
#' inp <- setInfasta(inp, aa)
#' 
#' # run signalp2 on the initial file:
#' step1_sp2 <- signalp(inp, version = 2, 'euk', run_mode = "starter", paths = my_pa)
#' 
#' # run wolfpsort on the signalp output:
#' w <- wolfpsort(step1_sp2, 'fungi', my_pa)

wolfpsort <- function(input_obj, organism, paths){
  # check that inputs are valid
  if (is(input_obj, "CBSResult")) {} else {stop('input_object does not belong to CBSResult superclass')}
  if (length(getOutfasta(input_obj)) == 0) {stop('the input object contains empty out_fasta slot')}
  allowed_organisms <- c('plant', 'animal', 'fungi')
  if (!(organism %in% allowed_organisms)) {stop('input organism is not allowed or does not exist')}
  
  message("running WoLF PSORT locally...")
  
  fasta <- getOutfasta(input_obj)
  message(paste("Number of submitted sequences...", length(fasta)))

  out_tmp <- tempfile()
  Biostrings::writeXStringSet(fasta, out_tmp)
  full_pa <- as.character(paths %>% filter(tool == 'wolfpsort') %>% select(path))
  wolf <- system(paste(full_pa,  organism, '<', out_tmp), intern = TRUE)
  
  #parse wolf output
  clean_strings <- function(x, field){unlist((stringr::str_split(x, " ")))[c(field)]}
  
  gene_id <- sapply(X = wolf, field = 1, clean_strings, USE.NAMES = FALSE)
  localization <- sapply(X = wolf, field = 2, clean_strings, USE.NAMES = FALSE)
  
  #assemble result tibble with gene id and most probable sibsellular localisation
  wolf_tbl <- as_tibble(data.frame(gene_id, localization)) %>% filter(gene_id != '#') %>% filter(localization == 'extr')
  
  message(paste('Number of candidate sequences with extracellular localisation...', nrow(wolf_tbl)))

  #assemble wolf result object:
  out_obj <- WolfResult(in_fasta = fasta,
                        out_fasta = fasta[wolf_tbl$gene_id], #place holder
                        wolf_tibble = wolf_tbl)
  
  if (validObject(out_obj)) {return(out_obj)}
}


