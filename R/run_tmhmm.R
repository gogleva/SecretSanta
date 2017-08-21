#' tmhmm function
#'
#' This function calls local TMHMM, expects CBSresult class objects with populated mature_fasta slot.
#' To generate it run signalp first.
#' @param input_obj input object, an instance of SignalpResult class, \cr
#'                  input should contain mature_fasta; in the full-length proteins \cr
#'                  N-terminal signal peptide could be erroneously \cr
#'                  predicted as TM domain, avoid this
#' @param paths tibble with paths to external dependencies, generated with \code{\link{manage_paths}} function
#' @param TM  allowed number of TM domains in mature peptides, recommended value <= 1; use 0 for strict filtering             
#' @export
#' @examples 
#'           
#' my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
#' inp <- SignalpResult()
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), use.names = TRUE)
#' inp <- setInfasta(inp, aa)
#' s1_sp2 <- signalp(inp, version = 2, 'euk', run_mode = "starter", paths = my_pa)
#' tm <- tmhmm(s1_sp2, paths = my_pa, TM = 1)

tmhmm <- function(input_obj, paths, TM) {

  if (TM >= 2) {warning('Recommended TM threshold values for mature peprides is 1')}  
  # check that input object belongs to a valid class
  if (is(input_obj, "SignalpResult")) {} else {stop('input_object does not belong to SignalpResult class')}
  
  # check that input object contains non-empty mature fasta slot
  s <- getSlots(class(input_obj))

  if ('mature_fasta' %in% names(s)) {
    if (length(getMatfasta(input_obj)) == 0) {
      stop('the input object contains an empty mature_fasta slot')
    }
  } else {
      stop('the input object does not contain mature_fasta slot')}
  
  
  #----- Run tmhmm
  message("running TMHMM locally...")
  
  fasta <- getMatfasta(input_obj) 
  out_tmp <- tempfile()
  Biostrings::writeXStringSet(fasta, out_tmp)
  
  message(paste('Submitted sequences...', length(fasta)))

  full_pa <- as.character(paths %>% dplyr::filter(tool == 'tmhmm') %>% dplyr::select(path))
  con <- system(paste(full_pa, out_tmp, '--short'), intern = TRUE)
  con_tmp <- tempfile()
  write(con, con_tmp)
  tm <- suppressMessages(readr::read_delim(con_tmp, '\t', col_names = F))
  
  names(tm) <- c("gene_id", "length", "ExpAA",
                     "First60", "PredHel", "Topology")
  
  # clean output values remove '... =' value
  clean_outp <- function(x) {unlist(stringr::str_split(x, '='))[2]}
  
  tm <- dplyr::mutate(tm,
                      length = sapply(tm$length, clean_outp, USE.NAMES = FALSE),
                      ExpAA = sapply(tm$ExpAA, clean_outp, USE.NAMES = FALSE),
                      First60 = sapply(tm$First60, clean_outp, USE.NAMES = FALSE),
                      PredHel = sapply(tm$PredHel, clean_outp, USE.NAMES = FALSE),
                      Topology = sapply(tm$Topology, clean_outp, USE.NAMES = FALSE)
                      )
  
  # change this lines in accordance with TM_thershold
  tm <- (tm %>% dplyr::filter(PredHel <= TM))
  
  message(paste('Candidate sequences with signal peptides and 0 TM domains in mature sequence...', nrow(tm)))
  
  # helper function: crop long names for AAStringSet object, return character vector
  crop_names <- function(x){unlist(stringr::str_split(x, " "))[1]}
  
  #generate cropped names for input fasta
  full_fasta <- getInfasta(input_obj)
  cropped_names <- unname(sapply(names(full_fasta), crop_names))
  
  #replace long names with cropped names
  names(full_fasta) <- cropped_names
  
  #get ids of candidate secreted proteins
  candidate_ids <- tm %>% dplyr::select(gene_id) %>% unlist(use.names = FALSE)
  out_fasta_tm <- full_fasta[candidate_ids]
  
  out_obj <- TMhmmResult(in_fasta = getOutfasta(input_obj), # original in fasta, full length proteins
                         out_fasta = out_fasta_tm, # out fasta, full length proteins
                         in_mature_fasta = fasta,
                         out_mature_fasta = fasta[candidate_ids],
                         tm_tibble = tm)
  
  # clean TMP files before exiting:
  
  junk <- dir(pattern = 'TMHMM*')
  file.remove(junk) 
  
  if (validObject(out_obj)) {return(out_obj)}
  }

