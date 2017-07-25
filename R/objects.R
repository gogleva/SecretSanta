#' An S4 class to represent intermediate and final outputs of the signalp prediction step
#' 
#' @slot in_fasta        initial fasta file, do we need to drag it along?
#' @slot out_fasta       output fasta with only positive candidates, i.e sequences predicted to be secreted at this step
#' @slot mature_fasta    fasta with mature sequences
#' @slot sp_version      version of signalp used to generate this object
#' @slot sp_tibble       Object of class tibble, contains ... columns:
#' \itemize{
#'   \item gene_id - unique id of the sequence
#'   \item Cmax - max raw cleavage site score (C-score)
#'   \item Cpos - amino acid position with max C-score
#'   \item Ymax - max combied cleavage site score (Y-score)
#'   \item Ypos - amino acid position with max Y-score
#'   \item Smax - max signal peptide score
#'   \item Spos - amino acid position with max S-score 
#'   \item D - a weighted average of the mean S and the max. Y scores.
#'             This is the score that is used to discriminate signal peptides from non-signal peptides.
#'   \item Smean - the average S-score of the possible signal peptide (from position 1 to the position immediately before the
#'             maximal Y-score)
#'   \item Prediction - final desision o whether the protein is secreted or not (Y/N)
#'   \item Dmaxcut - 
#'   \item Networks used - 
#'   }
#' @examples 
#' a <- SignalpResult()
#' setInfasta(a, bb)

SignalpResult <- setClass(
                    "SignalpResult",
                    slots = list(in_fasta = "AAStringSet", 
                                 out_fasta = "AAStringSet", 
                                 mature_fasta = "AAStringSet", 
                                 sp_version = "numeric",
                                 sp_tibble = "tbl_df"),
                 
                    prototype = list(
                                in_fasta = AAStringSet(),
                                out_fasta = AAStringSet(),
                                mature_fasta = AAStringSet(),
                                sp_version = 2,
                                sp_tibble = tibble()
                    ),
                    
                  # test if the data is consistent
                    validity = function(object)
                    {
                             if (length(object@in_fasta) < length(object@out_fasta)) {
                              return("Number of output sequences is grater than the number of input sequences.")
                             } else if (length(object@mature_fasta != length(object@out_fasta))) {
                              return("Numbers of sequences in output_fasta and mature_fasta do not match.") 
                             } else if (!(sp_version < 2 | sp_version >= 5)) {
                              return("signalp version is invalid") 
                             }
                             return(TRUE)
                    }  
                    )

# define getters ad setters for SignalpResult object
setGeneric(name = "setInfasta",
                def = function(theObject, in_fasta)
                {
                      standardGeneric("setInfasta")    
                }  
                )

setMethod(f = "setInfasta",
                signature = "SignalpResult",
                definition = function(theObject, in_fasta)
                {
                      theObject@in_fasta <- in_fasta
                      validObject(theObject)
                      return(theObject)
                }
                )


a <- SignalpResult()


#' An S4 class to represent intermediate and final outputs of the targetp prediction step
#' 
#' @slot in_fasta        initial fasta file, do we need to drag it along?
#' @slot out_fasta       output fasta with only positive candidates, i.e. not associated with cell organells
#' @slot tp_tibble       standard tibble with outputs obtained from targetp
#' 

TargetpResult <- setClass("TargetpResult",
                          slots = list(in_fasta = "file",
                                       out_fasta = "file",
                                       tp_tibble = "file",
                                       sp_tibble = "tbl_df")
)


#' An S4 class to represent intermediate and final outputs of the TMHMM prediction step
#' 
#' @slot in_fasta         initial fasta file
#' @slot in_mature_fasta  input mature fasta, 
#' @slot out_fasta        output fasta with only positive full-length candidates, candidates without TM in mature sequences
#' @slot out_mature_fasta output mature, conatins mature sequences without TM domains
#' @slot tm_tibble        tibble with outputs obtained from TMHMM
#' 

TMhmmResult <- setClass("TMhmmResult",
                           slots = list(in_fasta = "file",
                                        in_mature_fasta = "file",
                                        out_fasta = "file",
                                        tp_tibble = "file",
                                        sp_tibble = "tbl_df")
)


#' An S4 class to represent intermediate and final outputs of the ER motif checks (KDEL/HDEL)
#' 
#' @slot in_fasta         initial fasta file
#' @slot out_fasta        output fasta without sequences with terminal ER-retention signals
#' @slot er_list          list with seq ids containing terminal ER-retention motif
#' 

ErResult <- setClass("ErResult",
                        slots = list(in_fasta = "file",
                                     out_fasta = "file",
                                     er_list = "list")
)
