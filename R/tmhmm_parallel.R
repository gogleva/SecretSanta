# HELPER FUNCTIONS FOR tmhmm_parallel

#combine_TMhmmResult function
#'
#' This function combines multiple instances of TMhmmResult class, typically generated with parLapply
#' @param arguments - a list of TMhmmResult objects to be combined in one
#' @export
#' @examples 
#' inp2 <- CBSResult(in_fasta = readAAStringSet(system.file("extdata", "tail_prot.fasta", package = "SecretSanta")))
#' inp4 <- CBSResult(in_fasta = readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta")))
#' 
#' sp1 <- signalp(input_obj = inp2, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa)
#' sp3 <- signalp(input_obj = inp4, version = 2, organism_type = 'euk', run_mode = 'starter', paths = my_pa)
#' 
#' tm1 <- tmhmm(sp1, paths = my_pa, TM = 1)
#' tm3 <- tmhmm(sp3, paths = my_pa, TM = 1)

#' obj <- list(tm1, tm3)
#'  combined_tm <- combine_TMhmmResult(obj)


combine_TMhmmResult <- function(arguments) {
  if (all(sapply(arguments, is, 'TMhmmResult'))) {
  } else {                               
    stop('Some objects from arguments list do not belong to TMhmmResult class.')
  }
  
  c_in_fasta <- do.call(c, (lapply(arguments, getInfasta)))
  c_out_fasta <- do.call(c, (lapply(arguments, getOutfasta)))
  c_in_mature_fasta <- do.call(c, (lapply(arguments, getInMatfasta)))
  c_out_mature_fasta <- do.call(c, (lapply(arguments, getOutMatfasta)))
  c_tm_tibble <- do.call(rbind, (lapply(arguments, getTMtibble)))
  
  
  c_obj <- TMhmmResult(in_fasta = c_in_fasta,
                       out_fasta = c_out_fasta,
                       in_mature_fasta = c_in_mature_fasta,
                       out_mature_fasta = c_out_mature_fasta,
                       tm_tibble = c_tm_tibble)

}


# parallel version of TMHMM

#' tmhmm function
#' 
#' This function calls local TMHMM, expects CBSresult class objects with populated mature_fasta slot.
#' To generate CBSResult objects with mature_fasta slots run preffered version of signalp first.
#' Parallel version
#' @param input_obj input object, an instance of CBSResult class, \cr
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
#' s1_sp2_10K <- signalp(inp, version = 2, 'euk', run_mode = "starter", paths = my_pa)
#' tm <- tmhmm_parallel(s1_sp2, paths = my_pa, TM = 1)

tmhmm_parallel <- function(input_obj, paths, TM) {
  
  # helper function: crop long names for AAStringSet object, return character vector
  crop_names <- function(x){unlist(stringr::str_split(x, " "))[1]}
  
  # ----- Check that inputs are valid:
  
  if (TM >= 2) {warning('Recommended TM threshold values for mature peprides is 1')}  
  # check that input object belongs to a valid class
  if (is(input_obj, "CBSResult")) {} else {stop('input_object does not belong to CBSResult superclass')}
  
  # check that input object contains non-empty mature fasta slot
  s <- getSlots(class(input_obj))
  
  if ('mature_fasta' %in% names(s)) {
    if (length(getMatfasta(input_obj)) == 0) {
      stop('the input object contains an empty mature_fasta slot')
    }
  } else {
    stop('the input object does not contain mature_fasta slot')}
  
  
  # Crop names in the in_fasta filed of the input object:
  
  #generate cropped names for the input fasta ('full_fasta')
  full_fasta <- getInfasta(input_obj)
  cropped_names <- unname(sapply(names(full_fasta), crop_names))
  
  #replace long names with cropped names
  names(full_fasta) <- cropped_names
  
  simple_tmhmm <- function(fasta_tuple) {
  # fasta tuple - a lits of 2 elements: out_fasta, mature_fasta
  # (out_fasta, mature_fasta)
  
      
      #----- Run tmhmm
      message("running TMHMM locally...")
    
      outf <- fasta_tuple[[1]] #original out_fasta
      matf <- fasta_tuple[[2]] #original mature_fasta
     
      out_tmp <- tempfile()
      Biostrings::writeXStringSet(matf, out_tmp)
  
      message(paste('Number of submitted sequences...', length(matf)))
  
      full_pa <- as.character(paths %>% dplyr::filter(tool == 'tmhmm') %>% dplyr::select(path))
      tm <- tibble::as.tibble(read.table(text = (system(paste(full_pa, out_tmp, '--short'), intern = TRUE))))
      names(tm) <- c("gene_id", "length", "ExpAA",
                     "First60", "PredHel", "Topology")
  
      # helper function to clean output values remove '... =' value
      clean_outp <- function(x) {unlist(stringr::str_split(x, '='))[2]}
  
      tm <- dplyr::mutate(tm,
                          length = sapply(tm$length, clean_outp, USE.NAMES = FALSE),
                          ExpAA = sapply(tm$ExpAA, clean_outp, USE.NAMES = FALSE),
                          First60 = sapply(tm$First60, clean_outp, USE.NAMES = FALSE),
                          PredHel = sapply(tm$PredHel, clean_outp, USE.NAMES = FALSE),
                          Topology = sapply(tm$Topology, clean_outp, USE.NAMES = FALSE)
      )
  
     # select entries matching TM threshold:
     tm <- (tm %>% dplyr::filter(PredHel <= TM))
  
     message(paste('Number of candidate sequences with signal peptides and 0 TM domains in mature sequence...', nrow(tm)))
  
     #get ids of candidate secreted proteins
     candidate_ids <- tm %>% dplyr::select(gene_id) %>% unlist(use.names = FALSE)
     out_fasta_tm <- full_fasta[candidate_ids]
  
     out_obj <- TMhmmResult(in_fasta = outf, # original in fasta, full length proteins
                            out_fasta = out_fasta_tm, # out fasta, full length proteins
                            in_mature_fasta = matf,
                            out_mature_fasta = matf[candidate_ids],
                            tm_tibble = tm)
  
     # clean TMP files before exiting:
  
     junk <- dir(pattern = 'TMHMM*')
     file.remove(junk) 
     if (validObject(out_obj)) {return(out_obj)}
     
  }
  
  # estimate input file size
  
  all_mat_fasta <- getMatfasta(input_obj)
  all_out_fasta <- getOutfasta(input_obj)
  all_tup <- list(all_out_fasta, all_mat_fasta)
  
  if (length(all_mat_fasta) <= 10) {
    message('Ok for single file processing')
    return(simple_tmhmm(all_tup))
  } else {
    #split the files:
    
    split_out <- split_XStringSet(all_out_fasta, 10)
    split_mature <- split_XStringSet(all_mat_fasta, 10)
    fasta_tuples <- mapply(list, split_out, split_mature, SIMPLIFY=F)
    
    # do the parallel jobs
    
    # Calculate the number of cores
    no_cores <- detectCores()
    
    # Initiate cluster
    cl <- makeCluster(no_cores)
    # run parallel process
    
    clusterEvalQ(cl, library("SecretSanta"))
    clusterExport(cl=cl, varlist=c("my_pa")) # or path?
    
    result <- parLapply(cl, fasta_tuples, simple_tmhmm)
    stopCluster(cl)
    
    res_comb <- do.call(c,result)
    combined_TMhmmResult <- combine_TMhmmResult(unname(res_comb))
    
    tm_count <- nrow(getTMtibble(combined_TMhmmResult))    
    message(paste('Number of candidate sequences with less than', TM, 'TM domains...', tm_count))
    if (tm_count == 0) {warning(paste('TMHMM prediction yeilded 0 candidates with less than', TM, 'TM doamins'))}
    
    closeAllConnections()
    return(combined_TMhmmResult)
    #return(res_comb)
    }
}
