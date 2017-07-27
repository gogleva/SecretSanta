#' An S4 class to represent minimal basic structure of CBS output objects
#' 
#' @slot in_fasta        initial fasta file
#' @slot out_fasta       output fasta with only positive candidates, \cr
#'                       i.e those that passed tool filters

#' @examples 

CBSResult <- setClass("CBSResult",
                      
                slots = list(in_fasta = "AAStringSet",
                             out_fasta = "AAStringSet"),
                             

                prototype = list(in_fasta = AAStringSet(),
                                 out_fasta = AAStringSet()),
                
                validity = function(object)
                {
                  if (length(object@in_fasta) < length(object@out_fasta)) {
                    return("Number of output sequences is grater than the number of input sequences.")
                  } else if (any(grepl('[*$]', object@in_fasta))) {
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
#' @slot in_fasta        initial fasta file
#' @slot out_fasta       output fasta with only positive candidates, i.e sequences predicted to be secreted at this step
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
#'   \item D - a weighted average of the mean S and the max. Y scores.
#'             This is the score that is used to discriminate signal peptides from non-signal peptides.
#'   \item Smean - the average S-score of the possible signal peptide (from position 1 to the position immediately before #              the maximal Y-score)
#'   \item Prediction - final desision o whether the protein is secreted or not (Y/N)
#'   \item Dmaxcut - 
#'   \item Networks used - 
#'   }
#' @examples 
#' a <- SignalpResult()
#' setInfasta(a, bb)
#' a2 <- setInfasta(a, bb)
#' getInfasta(a2)
#  ss <- readAAStringSet("SecretSanta/inst/extdata/sample_prot_stop_codons.fasta", use.names = TRUE)
#  sc <- SignalpResult()

SignalpResult <- setClass(
                    "SignalpResult",
                    contains= "CBSResult",
                    slots = list(mature_fasta = "AAStringSet", 
                                 sp_version = "numeric",
                                 sp_tibble = "tbl_df"),
                 
                    prototype = list(mature_fasta = AAStringSet(),
                                sp_version = 2,
                                sp_tibble = tibble()
                    ),
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

#' An S4 class to represent intermediate and final outputs of the targetp prediction step
#' 
#' @slot tp_tibble       standard tibble with outputs obtained from targetp
#' 

TargetpResult <- setClass("TargetpResult",
                          contains = "CBSResult",
                          slots = list(tp_tibble = "tbl_df")
)


#' An S4 class to represent intermediate and final outputs of the TMHMM prediction step
#' 
#' @slot in_mature_fasta  input mature fasta, 
#' @slot out_mature_fasta output mature, conatins mature sequences without TM domains
#' @slot tm_tibble        tibble with outputs obtained from TMHMM
#' 

TMhmmResult <- setClass("TMhmmResult",
                           contains = "CBSResult",
                           slots = list(in_mature_fasta = "AAStringSet",
                                        out_mature_fasta = "AAStringSet",
                                        tm_tibble = "tbl_df")
)


#' An S4 class to represent intermediate and final outputs of the ER motif checks (KDEL/HDEL)
#' 
#' @slot er_list          list with seq ids containing terminal ER-retention motif
#' 

ErResult <- setClass("ErResult",
                        contains = "CBSResult",
                        slots = list(er_list = "list")
)
