#' An S4 class to represent intermediate and final outputs of the
#' secretome prediction pipeline
#' 
#' @slot paths           paths to CBS exceutibles, do I need it here?
#' @slot fasta           initial fasta file
#' @slot mature_fasta    fasta with mature sequences
#' @slot sp_tibble       standard tibble with outputs obtained from signalp
#' @slot tmhmm_tibble    tmhmm results
#' @slot targetp_tibble  targetp results
#' @slot er_motifs       ER motifs: (K/H)DEL
#' 
#' 
resultCookie <- setClass("resultCookie",
                    slots = list(paths = 'file',
                                 fasta = "file",
                                 mature_fasta = "file",
                                 sp_tibble = 'tbl_df',
                                 tmhmm_tibble = 'tbl_df',
                                 tagetp_tibble = 'tbl_df',
                                 er_motifs = 'character')
                        )




