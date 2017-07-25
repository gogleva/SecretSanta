#' An S4 class to represent intermediate and final outputs of the signalp prediction step
#' 
#' @slot paths           paths to CBS exceutibles, do I need it here?
#' @slot in_fasta        initial fasta file, do we need to drag it along?
#' @slot out_fasta       output fasta with only positive candidates, i.e sequences predicted to be secreted at this step
#' @slot mature_fasta    fasta with mature sequences
#' @slot sp_version      version of signalp used to generate this object
#' @slot sp_tibble       standard tibble with outputs obtained from signalp
#' 

SignalpResult <- setClass("SignalpResult",
                    slots = list(paths = 'file',
                                 in_fasta = "file",
                                 out_fasta = "file",
                                 mature_fasta = "file",
                                 sp_tibble = 'tbl_df')
                        )

#' An S4 class to represent intermediate and final outputs of the targetp prediction step
#' 
#' @slot paths           paths to CBS exceutibles, do I need it here?
#' @slot in_fasta        initial fasta file, do we need to drag it along?
#' @slot out_fasta       output fasta with only positive candidates, i.e. not associated with cell organells
#' @slot tp_tibble       standard tibble with outputs obtained from signalp
#' 

TargetppResult <- setClass("TargetpResult",
                          slots = list(paths = 'file',
                                       in_fasta = "file",
                                       out_fasta = "file",
                                       tp_tibble = "file",
                                       sp_tibble = 'tbl_df')
)




                    

