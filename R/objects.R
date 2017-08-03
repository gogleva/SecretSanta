#' An S4 class to represent minimal basic structure of CBS output objects.
#' 
#' @slot in_fasta        initial fasta file
#' @slot out_fasta       output fasta with only positive candidates,
#'                       i.e those that passed tool filters
#' @export CBSResult
#' @examples
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), use.names = TRUE)
#' 
#' # create an emty instance of CBSResult class
#' cbs <- CBSResult()
#' 
#' # populate in_fasta attribute with aa
#' cbs <- setInfasta(cbs, aa)
#' 
#' # populate out_fasta attribute with aa
#' cbs <- setOutfasta(cbs, aa)
#' 
#' # check that created CBSResult instance is valid
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
                  if (length(object@in_fasta) < length(object@out_fasta)) {
                    return("Number of output sequences is grater than the number of input sequences.")
                  }
                  
                  if (any(grepl('[*$]', object@in_fasta))) {
                    return("Input fasta contains stop codon symbols '*', please remove them.") 
                  }
                  
                  return(TRUE)
                }  
)
                                 

# Define accessors for CBSResult objects

# setter for in_fasta

setGeneric(name = "setInfasta",
           def = function(theObject, in_fasta)
           {
             standardGeneric("setInfasta")    
           }  
)


setMethod(f = "setInfasta",
          signature = "CBSResult",
          definition = function(theObject, in_fasta)
          {
            theObject@in_fasta <- in_fasta
            validObject(theObject)
            return(theObject)
          }
)

# getter for in_fasta

setGeneric(name = "getInfasta",
           def = function(theObject)
           {
             standardGeneric("getInfasta")    
           }  
)

setMethod(f = "getInfasta",
          signature = "CBSResult",
          definition = function(theObject)
          {
            return(theObject@in_fasta)
          }
)

# setter for out_fasta
setGeneric(name = "setOutfasta",
           def = function(theObject, out_fasta)
           {
             standardGeneric("setOutfasta")    
           }  
)

setMethod(f = "setOutfasta",
          signature = "CBSResult",
          definition = function(theObject, out_fasta)
          {
            theObject@out_fasta <- out_fasta
            validObject(theObject)
            return(theObject)
          }
)
# getter for out_fasta
setGeneric(name = "getOutfasta",
           def = function(theObject)
           {
             standardGeneric("getOutfasta")    
           }  
)

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
#' @examples 
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta", package = "SecretSanta"), use.names = TRUE)
#' 
#' # create an emty instance of CBSResult class
#' sr <- SignalpResult()
#' 
#' # populate in_fasta attribute with aa
#' sr <- setInfasta(sr, aa)
#' 
#' # run signalpeptide prediction and created fully pipulated instance of SignalpResult class
#' my_pa <- manage_paths(system.file("extdata", "sample_paths", package = "SecretSanta"))
#' step1_sp2 <- signalp(sr, version = 2, 'euk', run_mode = "starter", paths = my_pa)
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
                                sp_version = 2,
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
setGeneric(name = "setMatfasta",
           def = function(theObject, mature_fasta)
           {
             standardGeneric("setMatfasta")    
           }  
)

setMethod(f = "setMatfasta",
          signature = "SignalpResult",
          definition = function(theObject, mature_fasta)
          {
            theObject@mature_fasta <- mature_fasta
            validObject(theObject)
            return(theObject)
          }
)
# getter for mature_fasta
setGeneric(name = "getMatfasta",
           def = function(theObject)
           {
             standardGeneric("getMatfasta")    
           }  
)

setMethod(f = "getMatfasta",
          signature = "SignalpResult",
          definition = function(theObject)
          {
            return(theObject@mature_fasta)
          }
)

# setter for sp_version
setGeneric(name = "setSPversion",
           def = function(theObject, sp_version)
           {
             standardGeneric("setSPversion")    
           }  
)

setMethod(f = "setSPversion",
          signature = "SignalpResult",
          definition = function(theObject, sp_version)
          {
            theObject@sp_version <- sp_version
            validObject(theObject)
            return(theObject)
          }
)
# getter for sp_version
setGeneric(name = "getSPversion",
           def = function(theObject)
           {
             standardGeneric("getSPversion")    
           }  
)

setMethod(f = "getSPversion",
          signature = "SignalpResult",
          definition = function(theObject)
          {
            return(theObject@sp_version)
          }
)
# setter for sp_tible
setGeneric(name = "setSPtibble",
           def = function(theObject, sp_tibble)
           {
             standardGeneric("setSPtibble")    
           }  
)

setMethod(f = "setSPtibble",
          signature = "SignalpResult",
          definition = function(theObject, sp_tibble)
          {
            theObject@sp_tibble <- sp_tibble
            validObject(theObject)
            return(theObject)
          }
)
# getter for sp_version
setGeneric(name = "getSPtibble",
           def = function(theObject)
           {
             standardGeneric("getSPtibble")    
           }  
)

setMethod(f = "getSPtibble",
          signature = "SignalpResult",
          definition = function(theObject)
          {
            return(theObject@sp_tibble)
          }
)

#' An S4 class to represent outputs of WolfPsort
#' 
#' @slot wolf_tibble       standard tibble with outputs obtained from wolfpsort
#' @export WolfResult

WolfResult <- setClass("WolfResult",
                          contains = "CBSResult",
                          slots = list(wolf_tibble = "tbl_df")
)


#' An S4 class to represent intermediate and final outputs of the TMHMM prediction step
#' 
#' @slot in_fasta original in fasta, full length proteins
#' @slot in_mature_fasta  input mature fasta, extracted from the input SignalpResult object
#' @slot out_fasta output fasta, full-length sequences without TM after signal peptide
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
                                        tm_tibble = "tbl_df")
)


# setter for tm_tible
setGeneric(name = "setTMtibble",
           def = function(theObject, tm_tibble)
           {
             standardGeneric("setTMtibble")    
           }  
)

setMethod(f = "setTMtibble",
          signature = "TMhmmResult",
          definition = function(theObject, tm_tibble)
          {
            theObject@tm_tibble <- tm_tibble
            validObject(theObject)
            return(theObject)
          }
)
# getter for tm_tibble
setGeneric(name = "getTMtibble",
           def = function(theObject)
           {
             standardGeneric("getTMtibble")    
           }  
)

setMethod(f = "getTMtibble",
          signature = "TMhmmResult",
          definition = function(theObject)
          {
            return(theObject@tm_tibble)
          }
)


# getter for in_mature_fasta
setGeneric(name = "getInMatfasta",
           def = function(theObject)
           {
             standardGeneric("getInMatfasta")    
           }  
)

setMethod(f = "getInMatfasta",
          signature = "TMhmmResult",
          definition = function(theObject)
          {
            return(theObject@in_mature_fasta)
          }
)

# getter for out_mature_fasta
setGeneric(name = "getOutMatfasta",
           def = function(theObject)
           {
             standardGeneric("getOutMatfasta")    
           }  
)

setMethod(f = "getOutMatfasta",
          signature = "TMhmmResult",
          definition = function(theObject)
          {
            return(theObject@out_mature_fasta)
          }
)

#' An S4 class to represent intermediate and final outputs of the TMHMM prediction step
#'
#' @slot retained
#' @export ErResult

ErResult <- setClass("ErResult",
                      contains = "CBSResult",
                      slots = list(retained = "AAStringSet")
                     )

