#' tmhmm function
#'
#' This function calls local TMHMM
#' @param input_obj input object, an instance of CBSResult class, \cr
#'                  input should contain mature_fasta; in the full-length proteins \cr
#'                  N-terminal signal peptide could be erroneously \cr
#'                  predicted as TM domain, avoid this
#' @param paths tibble with paths to external dependencies, generated with \code{\link{manage_paths}} function                
#' @export
#' @examples 
#' 

tmhmm <- function(input_obj, paths) {
  
  # check that input object belongs to a valid class
  s <- getSlots(class(input_obj))

  if ('mature_fasta' %in% names(s)) {
    if (length(getMatfasta(input_obj)) == 0) {
      stop('the input object contains empty mature_fasta slot')
    }
  } else {
      stop('the input object does not contain mature_fasta slot')}
  
  message("running TMHMM locally...")
  
  fasta <- getMatfasta(input_obj) 
  out_tmp <- tempfile()
  Biostrings::writeXStringSet(fasta, out_tmp)
  
  full_pa <- as.character(secret_paths %>% dplyr::filter(tool == 'tmhmm') %>% dplyr::select(path))
  tm <- tibble::as.tibble(read.table(text = (system(paste(full_pa, out_tmp, '--short'), intern = TRUE))))
  names(tm) <- c("gene_id", "length", "ExpAA",
                     "First60", "PredHel", "Topology")
  tm <- (tm %>% dplyr::filter(PredHel == 'PredHel=0'))
  
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
  
  if (validObject(out_obj)) {return(out_obj)}
  }


### tests

# t2 <- tmhmm(step1_sp2) # obj of SignalpResult class
# t3 <- tmhmm(aa)
# t4 <- tmhmm(inp) 
# t22 <- tmhmm(step2_sp3)


