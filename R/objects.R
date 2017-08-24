#' An S4 class to represent minimal basic structure of CBS output objects.
#' 
#' @slot in_fasta        initial fasta file
#' @slot out_fasta       output fasta with only positive candidates,
#'                       i.e those that passed tool filters
#' @export CBSResult
#' @rdname CBS_methods
#' @examples
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata", 
#'                                   "sample_prot_100.fasta",
#'                                    package = "SecretSanta"),
#'                       use.names = TRUE)
#' 
#' # create an emty instance of CBSResult class
#' cbs <- CBSResult()
#' 
#' # populate in_fasta and out_fasta slots with aa
#' cbs <- CBSResult(in_fasta = aa, out_fasta = aa)
#' 
#' # populate out_fasta attribute with aa
#' cbs <- setOutfasta(cbs, aa)
#' 
#' # check that CBSResult instance is valid
#' validObject(cbs)
#' 
#' # extract in_fasta attribute
#' getInfasta(cbs)
#' 
#' # extract out_fasta attribute
#' getOutfasta(cbs)

CBSResult <- setClass("CBSResult",
                      
                slots = list(in_fasta = "AAStringSet",
                             out_fasta = "AAStringSet"),

                prototype = list(in_fasta = Biostrings::AAStringSet(),
                                 out_fasta = Biostrings::AAStringSet()),
                
                validity = function(object)
                {
                  #check that input is not dna or rna
                  
                  al <- Biostrings::alphabetFrequency(object@in_fasta)
                  
                                
                  if (nrow(al) == 1) {
                       most_frequent <- names(rev(sort(al[, colSums(al != 0) > 0])))[1:4] 
                   } else {
                       most_frequent <- names(rev(sort(colSums(al[, colSums(al != 0) > 0])))[1:4])
                   }
                  # 
                  
                  if (all(c("A", "C", "G", "T") %in% most_frequent)) {
                    return("Input sequence is DNA, please provide amino acid sequence.")
                  }
                  
                  if (all(c("A", "C", "G", "U") %in% most_frequent)) {
                    return("Input sequence is RNA, please provide amino acid sequence.")
                  }
                  
                  
                  if (length(object@in_fasta) < length(object@out_fasta)) {
                    return("Number of output sequences is grater than the number of input sequences.")
                  }
                  
                  if (any(grepl('[*$]', object@in_fasta))) {
                    return("Input fasta contains stop cteodon symbols '*', please remove them.") 
                  }
                  
                  if (any(duplicated(names(object@in_fasta)))) {
                    return("Duplicated gene ids in in_fasta")
                  }
                  
                  if (any(duplicated(names(object@out_fasta)))) {
                    return("Duplicated gene ids in out_fasta")
                  }
                  
                  return(TRUE)
                }  
)
                                 

# Define accessors for CBSResult objects

#' Accessors for CBSResult objects
#' @param theObject CBSResult object
#' @param in_fasta input fasta, AAStringSet object
#' @param out_fasta output fasta, AAStringSet object
#' @export
#' @docType methods
#' @rdname CBS_methods

setGeneric(name = "setInfasta",
           def = function(theObject, in_fasta)
           {
             standardGeneric("setInfasta")    
           }  
)

#' @export
#' @rdname  CBS_methods
#' @aliases setInfasta

setMethod(f = "setInfasta",
          signature = "CBSResult",
          definition = function(theObject, in_fasta)
          {
            theObject@in_fasta <- in_fasta
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "getInfasta",
           def = function(theObject)
           {
             standardGeneric("getInfasta")    
           }  
)

#' @export
#' @rdname  CBS_methods
#' @aliases getInfasta

setMethod(f = "getInfasta",
          signature = "CBSResult",
          definition = function(theObject)
          {
            return(theObject@in_fasta)
          }
)

setGeneric(name = "setOutfasta",
           def = function(theObject, out_fasta)
           {
             standardGeneric("setOutfasta")    
           }  
)

#' @export
#' @rdname  CBS_methods
#' @aliases setOutfasta

setMethod(f = "setOutfasta",
          signature = "CBSResult",
          definition = function(theObject, out_fasta)
          {
            theObject@out_fasta <- out_fasta
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "getOutfasta",
           def = function(theObject)
           {
             standardGeneric("getOutfasta")    
           }  
)

#' @export
#' @rdname  CBS_methods
#' @aliases getOutfasta

setMethod(f = "getOutfasta",
          signature = "CBSResult",
          definition = function(theObject)
          {
            return(theObject@out_fasta)
          }
)

#' An S4 class to represent intermediate and final outputs of the signalp prediction step
#' 
#' @slot mature_fasta    fasta with mature sequences
#' @slot sp_version      version of signalp used to generate this object
#' @slot sp_tibble       Object of class tibble, contains ... columns:
#' \itemize{
#'   \item gene_id - unique id of the sequence
#'   \item Cmax - max raw cleavage site score (C-score)
#'   \item Cpos - amino acid position with max C-score
#'   \item Cparsed - sp2 and sp3: remove this field?
#'   \item Ymax - max combied cleavage site score (Y-score)
#'   \item Ypos - amino acid position with max Y-score
#'   \item Smax - max signal peptide score
#'   \item Spos - amino acid position with max S-score 
#'   \item Smean - the average S-score of the possible signal peptide (from position 1 to the position immediately before #              the maximal Y-score)
#'   \item Prediction - final desision on whether the protein is secreted or not 
#'   }
#' @export SignalpResult
#' @rdname SignalpResult_methods 
#' @examples 
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata",
#'                                   "sample_prot_100.fasta",
#'                                   package = "SecretSanta"),
#'                       use.names = TRUE)
#' 
#' # create an emty instance of CBSResult class
#' sr <- SignalpResult()
#' 
#' # populate in_fasta attribute with aa
#' sr <- setInfasta(sr, aa)
#' 
#' # run signalpeptide prediction and create fully populated instance of SignalpResult class
#' my_pa <- manage_paths(system.file(
#'                                   "extdata",
#'                                   "sample_paths",
#'                                    package = "SecretSanta")
#'                                    )
#' step1_sp2 <- signalp(sr,
#'                      version = 2,
#'                      'euk',
#'                      run_mode = "starter",
#'                      paths = my_pa)
#' 
#' # access specific slots:
#' getOutfasta(step1_sp2)
#' getInfasta(step1_sp2)
#' getSPtibble(step1_sp2)
#' getSPversion(step1_sp2)
#' getMatfasta(step1_sp2) 


SignalpResult <- setClass(
                    "SignalpResult",
                    contains= "CBSResult",
                    slots = list(mature_fasta = "AAStringSet", 
                                 sp_version = "numeric",
                                 sp_tibble = "tbl_df"
                                 ),
                 
                    prototype = list(mature_fasta = Biostrings::AAStringSet(),
                                sp_version = numeric(0),
                                sp_tibble = tibble::tibble()
                                ),
                    
                    validity = function(object)
                    {
                      # check that mature sequences are shorter than full-length sequences
                      if (sum(width(object@out_fasta)) < sum(width(object@mature_fasta))) {
                        return("Mature sequences can not be shorter that full length sequences")
                      }
                      
                      # check that there are no duplicated gene ids in sp_tibble
                      if (nrow(object@sp_tibble) > 0) {
                        if (any(duplicated(object@sp_tibble$gene_id))) {
                        return("Duplicated gene ids in sp_tibble! ")}
                      }
                      
                      # check that all ids of mature_fasta are identical to ids in out_fasta
                      
                      if (!(identical(names(object@mature_fasta), names(object@out_fasta)))) {
                         return("Out_fasta ids do not match mature_fasta ids")
                        }
                      
                      # check that ids of mature_fasta are present in in_fasta
                      if (!(all(names(object@mature_fasta) %in% names(object@in_fasta)))) {
                        return("Out_fasta ids do not match in_fasta ids")
                      }
                      
                      # check that there are no zero length mature peptides - or make it a warning?
                      
                     if (any(width(object@mature_fasta) == 0)) {
                       return('mature fasta contains sequences of 0 length')
                     }
                    } 
                    )

# define accessor functions for SignalpResult object
# setter for mature_fasta

#' Accessors for SignalpResult objects
#' @param theObject SignalpResult object
#' @param mature_fasta seqeunces with clipped signal peptides, AAStringSet object
#' @param sp_version version of signalp used to generate SignalpResult object
#' @param sp_tibble parsed signalp output in tabular format
#' @export
#' @docType methods
#' @rdname SignalpResult_methods

setGeneric(name = "setMatfasta",
           def = function(theObject, mature_fasta)
           {
             standardGeneric("setMatfasta")    
           }  
)

#' @export
#' @rdname  SignalpResult_methods
#' @aliases setMatfasta

setMethod(f = "setMatfasta",
          signature = "SignalpResult",
          definition = function(theObject, mature_fasta)
          {
            theObject@mature_fasta <- mature_fasta
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "getMatfasta",
           def = function(theObject)
           {
             standardGeneric("getMatfasta")    
           }  
)

#' @export
#' @rdname  SignalpResult_methods
#' @aliases getMatfasta

setMethod(f = "getMatfasta",
          signature = "SignalpResult",
          definition = function(theObject)
          {
            return(theObject@mature_fasta)
          }
)

setGeneric(name = "setSPversion",
           def = function(theObject, sp_version)
           {
             standardGeneric("setSPversion")    
           }  
)

#' @export
#' @rdname  SignalpResult_methods
#' @aliases setSPversion

setMethod(f = "setSPversion",
          signature = "SignalpResult",
          definition = function(theObject, sp_version)
          {
            theObject@sp_version <- sp_version
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "getSPversion",
           def = function(theObject)
           {
             standardGeneric("getSPversion")    
           }  
)

#' @export
#' @rdname  SignalpResult_methods
#' @aliases getSPversion

setMethod(f = "getSPversion",
          signature = "SignalpResult",
          definition = function(theObject)
          {
            return(theObject@sp_version)
          }
)

setGeneric(name = "setSPtibble",
           def = function(theObject, sp_tibble)
           {
             standardGeneric("setSPtibble")    
           }  
)

#' @export
#' @rdname  SignalpResult_methods
#' @aliases setSPtibble

setMethod(f = "setSPtibble",
          signature = "SignalpResult",
          definition = function(theObject, sp_tibble)
          {
            theObject@sp_tibble <- sp_tibble
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "getSPtibble",
           def = function(theObject)
           {
             standardGeneric("getSPtibble")    
           }  
)

#' @export
#' @rdname  SignalpResult_methods
#' @aliases getSPtibble

setMethod(f = "getSPtibble",
          signature = "SignalpResult",
          definition = function(theObject)
          {
            return(theObject@sp_tibble)
          }
)

#' An S4 class to represent outputs of WolfPsort
#' 
#' @slot wolf_tibble       tibble with outputs obtained from wolfpsort
#' \itemize{
#'   \item gene_id - unique sequence id
#'   \item localization - the most probable localization
#' }
#' @export WolfResult

WolfResult <- setClass("WolfResult",
                          contains = "CBSResult",
                          slots = list(wolf_tibble = "tbl_df")
)

#' Accessors for WolfResult objects
#' @param theObject an object of WolfResult class
#' @param wolf_tibble parsed output of wolfpsort in tabular format
#' @export
#' @docType methods
#' @rdname WolfResult_methods

setGeneric(name = "setWOLFtibble",
           def = function(theObject, wolf_tibble)
           {
             standardGeneric("setWOLFtibble")    
           }  
)

#' @export
#' @rdname  WolfResult_methods
#' @aliases setWOLFtibble

setMethod(f = "setWOLFtibble",
          signature = "WolfResult",
          definition = function(theObject, wolf_tibble)
          {
            theObject@wolf_tibble <- wolf_tibble
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "getWOLFtibble",
           def = function(theObject)
           {
             standardGeneric("getWOLFtibble")    
           }  
)

#' @export
#' @rdname  WolfResult_methods
#' @aliases getWOLFtibble

setMethod(f = "getWOLFtibble",
          signature = "WolfResult",
          definition = function(theObject)
          {
            return(theObject@wolf_tibble)
          }
)

#' An S4 class to represent intermediate and final outputs of the TMHMM prediction step
#' 
#' @slot in_mature_fasta  input mature fasta, extracted from the input SignalpResult object
#' @slot out_mature_fasta output mature, conatins mature sequences without TM domains
#' @slot tm_tibble        tibble with outputs obtained from TMHMM
#' \itemize{
#'   \item gene_id - unique id of the sequence
#'   \item length - length of the protein sequence
#'   \item ExpAA - the expected number of amino acids intransmembrane helices
#'   \item First60 - the expected number of amino acids in transmembrane helices in the first 60 amino acids of the protein
#'   \item PredHel - the number of predicted transmembrane helices by N-best
#'   \item Topology - the topology predicted by N-best' 
#' }
#' @export TMhmmResult

TMhmmResult <- setClass("TMhmmResult",
                           contains = "CBSResult",
                           slots = list(in_mature_fasta = "AAStringSet",
                                        out_mature_fasta = "AAStringSet",
                                        tm_tibble = "tbl_df"),
                        
                           validity = function(object) {   
                          
                           # check that there are o duplicated ids in the input and output fastas and tm_tibble
                             if (nrow(object@tm_tibble) > 0) {
                                if (any(duplicated(object@tm_tibble$gene_id))) {
                                  return("Duplicated gene ids in sp_tibble! ")}
                             } 
                             
                             if (any(duplicated(names(object@in_mature_fasta)))) {
                                  return("Duplicated gene ids in in_mature_fasta")
                             }
                             
                             if (any(duplicated(names(object@out_mature_fasta)))) {
                               return("Duplicated gene ids in out_mature_fasta")
                             }
                             
                           # check that ids in out_mature_fasta match ids in out_fasta
                              
                               if (!(identical(names(object@out_mature_fasta), names(object@out_fasta)))) {
                                 return("out_fasta ids do not match out_mature_fasta ids")
                               }
                            }  
)


#' Accessors for TMhmmResul objects
#' @param theObject object of TMhmmREsult class
#' @param tm_tibble parsed tmhmm output in tabular format
#' @export
#' @docType methods
#' @rdname TMhmmResult_methods

setGeneric(name = "setTMtibble",
           def = function(theObject, tm_tibble)
           {
             standardGeneric("setTMtibble")    
           }  
)
#' @export
#' @rdname  TMhmmResult_methods
#' @aliases setTMtibble
setMethod(f = "setTMtibble",
          signature = "TMhmmResult",
          definition = function(theObject, tm_tibble)
          {
            theObject@tm_tibble <- tm_tibble
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "getTMtibble",
           def = function(theObject)
           {
             standardGeneric("getTMtibble")    
           }  
)
#' @export
#' @rdname  TMhmmResult_methods
#' @aliases getTMtibble
setMethod(f = "getTMtibble",
          signature = "TMhmmResult",
          definition = function(theObject)
          {
            return(theObject@tm_tibble)
          }
)

setGeneric(name = "getInMatfasta",
           def = function(theObject)
           {
             standardGeneric("getInMatfasta")    
           }  
)
#' @export
#' @rdname  TMhmmResult_methods
#' @aliases getInMatfasta
setMethod(f = "getInMatfasta",
          signature = "TMhmmResult",
          definition = function(theObject)
          {
            return(theObject@in_mature_fasta)
          }
)

setGeneric(name = "getOutMatfasta",
           def = function(theObject)
           {
             standardGeneric("getOutMatfasta")    
           }  
)
#' @export
#' @rdname  TMhmmResult_methods
#' @aliases getOutMatfasta
setMethod(f = "getOutMatfasta",
          signature = "TMhmmResult",
          definition = function(theObject)
          {
            return(theObject@out_mature_fasta)
          }
)

#' An S4 class to represent intermediate and final outputs of the TMHMM prediction step
#'
#' @slot retained - sequences with ER retention signals
#' @export ErResult

ErResult <- setClass("ErResult",
                      contains = "CBSResult",
                      slots = list(retained = "AAStringSet")
                     )


#' An S4 class to represent intermediate and final outputs of the targetP prediction step
#' 
#' @slot tp_tibble        tibble with outputs obtained from targetp
#' \itemize{
#' \item gene_id - unique id of the sequence
#' \item length - length of the protein sequence
#' \item mTP - mitochondrial NN score
#' \item sp - signal peptide NN score
#' \item other - any onther NN score
#' \item TP_localization - 	Prediction of localization, based on the scores; the possible values are:
#' } 
#'   \itemize{
#'      \item C	- Chloroplast, i.e. the sequence contains cTP, a chloroplast transit peptide;
#'      \item M	- Mitochondrion, i.e. the sequence contains mTP, a mitochondrial targeting peptide;
#'      \item S -	Secretory pathway, i.e. the sequence contains SP, a signal peptide;
#'      \item _ -	Any other location;
#'      \item "don't know" - indicates that cutoff restrictions were set (see instructions) and the winning network output score was below the requested cutoff for that category.
#'      \item RC - Reliability class, from 1 to 5, where 1 indicates the strongest prediction. RC is a measure of the size of the difference ('diff') between the highest (winning) and the second highest output scores. There are 5 reliability classes, defined as follows:
#'   \itemize{ 
#'     \item 1 - diff > 0.800;
#'     \item 2 - 0.800 > diff > 0.600
#'     \item 3 - 0.600 > diff > 0.400
#'     \item 4 - 0.400 > diff > 0.200
#'     \item 5 - 0.200 > diff
#'     }
#' }
#' @export TargetpResult

TargetpResult <- setClass("TargetpResult",
                        contains = "CBSResult",
                        slots = list(tp_tibble = "tbl_df"),
                        
                        validity = function(object) {   
                          
                          # check that there are o duplicated ids in the input and output fastas and tm_tibble
                          if (nrow(object@tp_tibble) > 0) {
                            if (any(duplicated(object@tp_tibble$gene_id))) {
                              return("Duplicated gene ids in sp_tibble! ")}
                          } 
                        }  
)

#' Accessors for TargetpResult objects
#' @param theObject an object of TargetpResult class
#' @param tp_tibble parsed targetp output in tabular format
#' @export
#' @docType methods
#' @rdname TargetpResult_methods

setGeneric(name = "setTPtibble",
           def = function(theObject, tp_tibble)
           {
             standardGeneric("setTPtibble")    
           }  
)

#' @export
#' @rdname  TargetpResult_methods
#' @aliases setTPtibble

setMethod(f = "setTPtibble",
          signature = "TargetpResult",
          definition = function(theObject, tp_tibble)
          {
            theObject@tp_tibble <- tp_tibble
            validObject(theObject)
            return(theObject)
          }
)

setGeneric(name = "getTPtibble",
           def = function(theObject)
           {
             standardGeneric("getTPtibble")    
           }  
)

#' @export
#' @rdname  TargetpResult_methods
#' @aliases getTPtibble

setMethod(f = "getTPtibble",
          signature = "TargetpResult",
          definition = function(theObject)
          {
            return(theObject@tp_tibble)
          }
)