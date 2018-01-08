#' accessor functions for objects of \code{\link{CBSResult}} S4 class, minimal parent class for
#' SecretSanat outputs
#'
#' @slot in_fasta        initial fasta file
#' @slot out_fasta       output fasta file with only positive candidates,
#'                       i.e those that passed tool filters
#' @export CBSResult
#' @return CBSResult object
#' @rdname CBS_methods
#' @examples
#' library(Biostrings)
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata","sample_prot_100.fasta",
#'                         package = "SecretSanta"))
#' # create an instance of CBSResult class:
#' cbs <- CBSResult( in_fasta = aa, out_fasta = aa)
#' # check that CBSResult instance is valid
#' validObject(cbs)
#'
#' # extract in_fasta attribute
#' getInfasta(cbs)
#'
#' # extract out_fasta attribute
#' getOutfasta(cbs)

CBSResult <- setClass("CBSResult",
    contains = "AAStringSetList",                  
    slots = c(in_fasta = "AAStringSet",
                 out_fasta = "AAStringSet",
                 seqList = "AAStringSetList")
    )

# define validity fuction for CBSResult class:

validCBSResult <- function(object) {
    al <- Biostrings::alphabetFrequency(object@seqList$in_fasta)
    
    if (nrow(al) == 1) {
        most_frequent <- names(rev(sort(al[ , colSums(al != 0) > 0])))[1:4]
    } else {
        most_frequent <-
            names(rev(sort(colSums(al[, colSums(al != 0) > 0])))[1:4])
    }
    
    if (all(c("A", "C", "G", "T") %in% most_frequent)) {
        return("Input sequence is DNA, please provide amino acid sequence.")
    }
    
    if (all(c("A", "C", "G", "U") %in% most_frequent)) {
        return("Input sequence is RNA, amino acid sequence required.")
    }
    
    
    if (length(object@seqList$in_fasta) < length(object@seqList$out_fasta)) {
        return("Number of output sequences > the number of input sequences.")
    }
    
    if (any(grepl('[*$]', object@seqList$in_fasta))) {
        return("Input fasta contains stop codon symbols '*'")
    }
    
    if (any(duplicated(names(object@seqList$in_fasta)))) {
        return("Duplicated gene ids in in_fasta")
    }
    
    if (any(duplicated(names(object@seqList$out_fasta)))) {
        return("Duplicated gene ids in out_fasta")
    }
    
    return(TRUE)
}
    
# set validity function
setValidity("CBSResult", validCBSResult)

# define custom constructor to pack in_fasta and out_fasta slots in 
# AAStringSetList

setMethod(f = 'initialize',
          signature = "CBSResult",
          definition = function(.Object,
                                in_fasta = AAStringSet(),
                                out_fasta = AAStringSet()) {
              .Object@seqList <- Biostrings::AAStringSetList('in_fasta' = in_fasta, 
                                                            'out_fasta' = out_fasta)
              validObject(.Object)
              .Object
          }
)

# define show method for CBSResult class:

setMethod("show",
          signature = 'CBSResult',
          definition = function(object){
              cat(paste("An object of class", class(object), "\n"))
              print(elementNROWS(object@seqList))
              
          })

# define accessors for CBSResult objects

#' accessors for CBSResult objects
#' @param theObject \code{\link{CBSResult}} object
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
            theObject@seqList$in_fasta <- in_fasta
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
            return(theObject@seqList$in_fasta)
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

setMethod(
    f = "setOutfasta",
    signature = "CBSResult",
    definition = function(theObject, out_fasta)
    {
        theObject@seqList$out_fasta <- out_fasta
        validObject(theObject)
        return(theObject)
    }
)

setGeneric(
    name = "getOutfasta",
    def = function(theObject)
    {
        standardGeneric("getOutfasta")
    }
)

#' @export
#' @rdname  CBS_methods
#' @aliases getOutfasta

setMethod(
    f = "getOutfasta",
    signature = "CBSResult",
    definition = function(theObject)
    {
        return(theObject@seqList$out_fasta)
    }
)

setGeneric(name = "getFastas",
           def = function(theObject)
           {
               standardGeneric("getFastas")
           }
)

#' @export
#' @rdname  CBS_methods
#' @aliases getFastas

setMethod(f = "getFastas",
          signature = "CBSResult",
          definition = function(theObject)
          {
              return(theObject@seqList)
          }
)

#' accessor functions for objects of SignalpResult S4 class, intermediate and
#' final outputs of the signalp prediction step
#'
#' @slot mature_fasta    fasta with mature sequences
#' @slot sp_version      version of \code{signalp} used to generate this object
#' @slot sp_tibble       Object of class tibble, columns:
#' \itemize{
#' \item gene_id - unique id of the sequence
#' \item Cmax - max raw cleavage site score (C-score)
#' \item Cpos - amino acid position with max C-score
#' \item Cparsed - sp2 and sp3: remove this field?
#' \item Ymax - max combied cleavage site score (Y-score)
#' \item Ypos - amino acid position with max Y-score
#' \item Smax - max signal peptide score
#' \item Spos - amino acid position with max S-score
#' \item Smean - the average S-score of the possible signal peptide
#' (from position 1 to the position immediately before the maximal Y-score)
#' \item Prediction - final desision on whether the protein is secreted
#' }
#' @export SignalpResult
#' @rdname SignalpResult_methods
#' @return SignalpResult object
#' @examples
#' # read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata","sample_prot_100.fasta",
#' package = "SecretSanta"))
#' # create an instance of CBSResult class
#' sr <- CBSResult(in_fasta = aa[1:20])
#' # run signalpeptide prediction and create fully populated instance of
#' # SignalpResult class
#' step1_sp2 <- signalp(sr, version = 2, organism = 'euk', run_mode = "starter")
#' # access specific slots:
#' getOutfasta(step1_sp2)
#' getInfasta(step1_sp2)
#' getSPtibble(step1_sp2)
#' getSPversion(step1_sp2)
#' getMatfasta(step1_sp2)


SignalpResult <- setClass(
    "SignalpResult",
    contains = "CBSResult",
    slots = list(
        mature_fasta = "AAStringSet",
        sp_version = "numeric",
        sp_tibble = "tbl_df"))

validSignalpResult <- function(object) {
    
        # check that mature sequences are shorter than full-length sequences
        if (sum(Biostrings::width(object@out_fasta)) <
            sum(Biostrings::width(object@mature_fasta))) {
            return("Mature sequences can not be shorter that full length ones")
        }

        # check that there are no duplicated gene ids in sp_tibble
        if (nrow(object@sp_tibble) > 0) {
            if (any(duplicated(object@sp_tibble$gene_id))) {
                return("Duplicated gene ids in sp_tibble! ")
            }
        }

        # check that all ids of mature_fasta are identical to ids
        # in out_fasta

        if (!(identical(names(object@mature_fasta),
                        names(object@out_fasta)))) {
            return("Out_fasta ids do not match mature_fasta ids")
        }
        

        # check that ids of mature_fasta are present in in_fasta
        if (!(all(names(object@mature_fasta) %in%
                    names(object@in_fasta)))) {
            return("Out_fasta ids do not match in_fasta ids")
        }
    
        # check that there are no zero length mature peptides

        if (any(width(object@mature_fasta) == 0)) {
            return('mature fasta contains sequences of 0 length')
        }
    }

# set validity function for SignalpResult objects
setValidity("SignalpResult", validSignalpResult)

# define constructor to add mature_fasta to seqList slot, together with in_fasta
# and out_fasta, oragnised AAStringSetList

setMethod(f = 'initialize',
          signature = "SignalpResult",
          definition = function(.Object,
                                mature_fasta = Biostrings::AAStringSet(),
                                sp_version = numeric(0),
                                sp_tibble = tibble::tibble()
                                ) {
              .Object@seqList <- Biostrings::AAStringSetList(
                  'in_fasta' = .Object@in_fasta,
                  'out_fasta' = .Object@out_fasta,
                  'mature_fasta' = mature_fasta)
              .Object@sp_version <- sp_version
              .Object@sp_tibble <- sp_tibble
              validObject(.Object)
              .Object
          }
         )

# define show method for SignalpResult class:

setMethod("show",
          signature = 'SignalpResult',
          definition = function(object){
              cat(paste("An object of class", class(object), "\n"))
              print(elementNROWS(object@seqList))
              cat(paste('SignalP version ... ', object@sp_version, "\n"))
              cat('SignalP tabular output:', '\n')
              print(object@sp_tibble)
              
          })

#' accessors for SignalpResult objects
#' @param theObject \code{\link{SignalpResult}} object
#' @param mature_fasta sequences with clipped signal peptides, AAStringSet
#' object
#' @param sp_version version of signalp used to generate SignalpResult object
#' @param sp_tibble parsed signalp output in tabular format
#' @export
#' @docType methods
#' @rdname SignalpResult_methods

setGeneric(
    name = "setMatfasta",
    def = function(theObject, mature_fasta)
    {
        standardGeneric("setMatfasta")
    }
)

#' @export
#' @rdname  SignalpResult_methods
#' @aliases setMatfasta

setMethod(
    f = "setMatfasta",
    signature = "SignalpResult",
    definition = function(theObject, mature_fasta)
    {
        theObject@seqList@mature_fasta <- mature_fasta
        validObject(theObject)
        return(theObject)
    }
)

setGeneric(
    name = "getMatfasta",
    def = function(theObject)
    {
        standardGeneric("getMatfasta")
    }
)

#' @export
#' @rdname  SignalpResult_methods
#' @aliases getMatfasta

setMethod(
    f = "getMatfasta",
    signature = "SignalpResult",
    definition = function(theObject)
    {
        return(theObject@seqList$mature_fasta)
    }
)

setGeneric(
    name = "setSPversion",
    def = function(theObject, sp_version)
    {
        standardGeneric("setSPversion")
    }
)

#' @export
#' @rdname  SignalpResult_methods
#' @aliases setSPversion

setMethod(
    f = "setSPversion",
    signature = "SignalpResult",
    definition = function(theObject, sp_version)
    {
        theObject@sp_version <- sp_version
        validObject(theObject)
        return(theObject)
    }
)

setGeneric(
    name = "getSPversion",
    def = function(theObject)
    {
        standardGeneric("getSPversion")
    }
)

#' @export
#' @rdname  SignalpResult_methods
#' @aliases getSPversion

setMethod(
    f = "getSPversion",
    signature = "SignalpResult",
    definition = function(theObject)
    {
        return(theObject@sp_version)
    }
)

setGeneric(
    name = "setSPtibble",
    def = function(theObject, sp_tibble)
    {
        standardGeneric("setSPtibble")
    }
)

#' @export
#' @rdname  SignalpResult_methods
#' @aliases setSPtibble

setMethod(
    f = "setSPtibble",
    signature = "SignalpResult",
    definition = function(theObject, sp_tibble)
    {
        theObject@sp_tibble <- sp_tibble
        validObject(theObject)
        return(theObject)
    }
)

setGeneric(
    name = "getSPtibble",
    def = function(theObject)
    {
        standardGeneric("getSPtibble")
    }
)

#' @export
#' @rdname  SignalpResult_methods
#' @aliases getSPtibble

setMethod(
    f = "getSPtibble",
    signature = "SignalpResult",
    definition = function(theObject)
    {
        return(theObject@sp_tibble)
    }
)

#' accessor functions for objects of WolfResult S4 class, outputs of the wolf
#' prediction step
#'
#' @slot wolf_tibble    tibble with outputs obtained from wolfpsort
#' \itemize{
#' \item gene_id - unique sequence id
#' \item localization - the most probable localization
#' }
#' @export WolfResult
#' @rdname  WolfResult_methods
#' @return WolfResult object
#' @examples
#' library(Biostrings)
#' aa <- readAAStringSet(system.file("extdata","sample_prot_100.fasta",
#' package = "SecretSanta"))
#' inp <- CBSResult(in_fasta = aa[1:10])
#' s1_sp2 <- signalp(inp, version = 2, organism = 'euk',
#' run_mode = "starter")
#' w <- wolfpsort(s1_sp2, 'fungi')
#' class(w)
#' #access result tibble:
#' getWOLFtibble(w)
#' #acess in_fasta and out_fasta slots:
#' getInfasta(w)
#' getOutfasta(w)

WolfResult <- setClass("WolfResult", 
                       contains = "CBSResult",
                       slots = list(wolf_tibble = "tbl_df"))

# constructor for WolfResult objects

setMethod(f = 'initialize',
          signature = "WolfResult",
          definition = function(.Object,
                                wolf_tibble = tibble::tibble()) {
              .Object@wolf_tibble <- wolf_tibble
          #     validObject(.Object)
              .Object
          }
)

# show method for WolfResult object
setMethod("show",
          signature = 'WolfResult',
          definition = function(object){
              cat(paste("An object of class", class(object), "\n"))
              if (length(object@seqList) != 0) {
              print(elementNROWS(object@seqList))
              } else {
              cat("all fasta slots are empty", '\n')          
                  }
              cat('WoLFPSort tabular output:', '\n')
              print(object@wolf_tibble)
              
          })

#' accessors for WolfResult objects
#' @param theObject an object of WolfResult class
#' @param wolf_tibble parsed output of wolfpsort in tabular format
#' @export
#' @docType methods
#' @rdname WolfResult_methods

setGeneric(
    name = "setWOLFtibble",
    def = function(theObject, wolf_tibble)
    {
        standardGeneric("setWOLFtibble")
    }
)

#' @export
#' @rdname  WolfResult_methods
#' @aliases setWOLFtibble

setMethod(
    f = "setWOLFtibble",
    signature = "WolfResult",
    definition = function(theObject, wolf_tibble)
    {
        theObject@wolf_tibble <- wolf_tibble
        validObject(theObject)
        return(theObject)
    }
)

setGeneric(
    name = "getWOLFtibble",
    def = function(theObject)
    {
        standardGeneric("getWOLFtibble")
    }
)

#' @export
#' @rdname  WolfResult_methods
#' @aliases getWOLFtibble

setMethod(
    f = "getWOLFtibble",
    signature = "WolfResult",
    definition = function(theObject)
    {
        return(theObject@wolf_tibble)
    }
)

#' accessor functions for objects of TMhmmResult S4 class, intermediate and
#' final outputs of the tmhmm prediction step
#' @slot in_mature_fasta  input mature fasta, extracted from the input
#' SignalpResult object
#' @slot out_mature_fasta output mature, conatins mature sequences without
#'  TM domains
#' @slot tm_tibble        tibble with outputs obtained from TMHMM
#' \itemize{
#' \item gene_id - unique id of the sequence
#' \item length - length of the protein sequence
#' \item ExpAA - the expected number of amino acids intransmembrane helices
#' \item First60 - the expected number of amino acids in transmembrane
#' helices in the first 60 amino acids of the protein
#' \item PredHel - the number of predicted transmembrane helices by N-best
#' \item Topology - the topology predicted by N-best'
#' }
#' @export TMhmmResult
#' @rdname TMhmmResult_methods
#' @examples
#' aa <- readAAStringSet(system.file("extdata", "sample_prot_100.fasta",
#' package = "SecretSanta"))
#' inp <- CBSResult(in_fasta = aa[1:10])
#' s1_sp2 <- signalp(inp, version = 2, organism = 'euk',
#' run_mode = "starter")
#' tm <- tmhmm(s1_sp2, TM = 1)
#' class(tm)
#' getTMtibble(tm)

TMhmmResult <- setClass(
    "TMhmmResult",
    contains = "CBSResult",
    slots = list(
        in_mature_fasta = "AAStringSet",
        out_mature_fasta = "AAStringSet",
        tm_tibble = "tbl_df"
    ),

    validity = function(object) {
        # check that there are o duplicated ids in the input
        # and output fastas and tm_tibble
        if (nrow(object@tm_tibble) > 0) {
            if (any(duplicated(object@tm_tibble$gene_id))) {
                return("Duplicated gene ids in sp_tibble! ")
            }
        }

        if (any(duplicated(names(object@in_mature_fasta)))) {
            return("Duplicated gene ids in in_mature_fasta")
        }

        if (any(duplicated(names(object@out_mature_fasta)))) {
            return("Duplicated gene ids in out_mature_fasta")
        }

        # check that ids in out_mature_fasta match ids in out_fasta

        if (!(identical(names(object@out_mature_fasta),
                        names(object@out_fasta)))) {
            return("out_fasta ids do not match out_mature_fasta ids")
        }
    }
)




#' Accessors for TMhmmResul objects
#' @param theObject object of TMhmmREsult class
#' @param tm_tibble parsed tmhmm output in tabular format
#' @export
#' @docType methods
#' @return TMhmmResult object
#' @rdname TMhmmResult_methods

setGeneric(
    name = "setTMtibble",
    def = function(theObject, tm_tibble)
    {
        standardGeneric("setTMtibble")
    }
)
#' @export
#' @rdname  TMhmmResult_methods
#' @aliases setTMtibble
setMethod(
    f = "setTMtibble",
    signature = "TMhmmResult",
    definition = function(theObject, tm_tibble)
    {
        theObject@tm_tibble <- tm_tibble
        validObject(theObject)
        return(theObject)
    }
)

setGeneric(
    name = "getTMtibble",
    def = function(theObject)
    {
        standardGeneric("getTMtibble")
    }
)
#' @export
#' @rdname  TMhmmResult_methods
#' @aliases getTMtibble
setMethod(
    f = "getTMtibble",
    signature = "TMhmmResult",
    definition = function(theObject)
    {
        return(theObject@tm_tibble)
    }
)

setGeneric(
    name = "getInMatfasta",
    def = function(theObject)
    {
        standardGeneric("getInMatfasta")
    }
)
#' @export
#' @rdname  TMhmmResult_methods
#' @aliases getInMatfasta
setMethod(
    f = "getInMatfasta",
    signature = "TMhmmResult",
    definition = function(theObject)
    {
        return(theObject@in_mature_fasta)
    }
)

setGeneric(
    name = "getOutMatfasta",
    def = function(theObject)
    {
        standardGeneric("getOutMatfasta")
    }
)
#' @export
#' @rdname  TMhmmResult_methods
#' @aliases getOutMatfasta
setMethod(
    f = "getOutMatfasta",
    signature = "TMhmmResult",
    definition = function(theObject)
    {
        return(theObject@out_mature_fasta)
    }
)

#' accessor functions for objects of ErResult S4 class, outputs of the
#' check_khedel function
#' @slot retained_fasta  sequences with ER retention signals
#' @export ErResult
#' @rdname ErResult_methods
#' @examples
#' aa <- readAAStringSet(system.file("extdata","small_prot.fasta",
#' package = "SecretSanta"), use.names = TRUE)
#' er <- check_khdel(CBSResult(in_fasta = aa), run_mode = 'starter')
#' class(er)

ErResult <- setClass("ErResult", contains = "CBSResult",
                        slots = list(retained_fasta = "AAStringSet"))


#' accessors for ErResult objects
#' @param theObject an object of ErResult class
#' @slot retained_fasta sequences with ER retention signals
#' @export
#' @docType methods
#' @rdname ErResult_methods

setGeneric(
    name = "getRetainedfasta",
    def = function(theObject)
    {
        standardGeneric("getRetainedfasta")
    }
)

#' @export
#' @rdname  ErResult_methods
#' @aliases getRetainedfasta


setMethod(
    f = "getRetainedfasta",
    signature = "ErResult",
    definition = function(theObject)
    {
        return(theObject@retained_fasta)
    }
)



#' accessor functions for objects of TargetpResult S4 class, intermediate and
#' final outputs of the targetp prediction step
#'
#' @slot tp_tibble        tibble with outputs obtained from targetp
#' \itemize{
#' \item    gene_id - unique id of the sequence
#' \item    length - length of the protein sequence
#' \item    mTP - mitochondrial NN score
#' \item    sp - signal peptide NN score
#' \item    other - any onther NN score
#' \item    TP_localization - Prediction of localization, based on the scores;
#' the possible values are:
#' }
#' \itemize{
#' \item    C - Chloroplast, i.e. the sequence contains cTP, a chloroplast
#' transit peptide;
#' \item    M - Mitochondrion, i.e. the sequence contains mTP, a mitochondrial
#' targeting peptide;
#' \item    S - Secretory pathway, i.e. the sequence contains SP, a signal
#' peptide;
#' \item    _ - Any other location;
#' \item    "don't know" - indicates that cutoff restrictions were set
#' (see instructions) and the winning network output score was below the
#' requested cutoff for that category.
#' \item    RC - Reliability class, from 1 to 5, where 1 indicates the
#' strongest prediction. RC is a measure of the size of the difference
#' ('diff') between #'  the highest (winning) and the second highest output
#' scores. There are 5 reliability classes, defined as follows:
#' \itemize{
#'     \item 1 - diff > 0.800;
#'     \item 2 - 0.800 > diff > 0.600
#'     \item 3 - 0.600 > diff > 0.400
#'     \item 4 - 0.400 > diff > 0.200
#'     \item 5 - 0.200 > diff
#'     }
#' }
#' @export TargetpResult
#' @rdname TargetpResult_methods
#' @examples
#' #read fasta file in AAStringSet object
#' aa <- readAAStringSet(system.file("extdata","sample_prot_100.fasta",
#' package = "SecretSanta"))
#' #assign this object to the input_fasta slot of SignalpResult object
#' inp <- CBSResult(in_fasta = aa[1:10])
#' tp_result <- targetp(input_obj = inp, network = 'N',
#' run_mode = 'starter')
#' class(tp_result)
#' getTPtibble(tp_result)
#' getInfasta(tp_result)
#' getOutfasta(tp_result)

TargetpResult <- setClass(
    "TargetpResult",
    contains = "CBSResult",
    slots = c(tp_tibble = "tbl_df"))

# validity function for TargetpResult class

validTargetpResult <- function(object) {
        # check that there are no duplicated ids in the input
        # and output fastas and tm_tibble
        if (nrow(object@tp_tibble) > 0) {
            if (any(duplicated(object@tp_tibble$gene_id))) {
                return("Duplicated gene ids in sp_tibble! ")
            }
        }
    }

# set validity function
setValidity("TargetpResult", validTargetpResult)

# constructor for TargetpResult objects:

setMethod(f = 'initialize',
          signature = "TargetpResult",
          definition = function(.Object,
                                tp_tibble = tibble::tibble()) {
              .Object@tp_tibble <- tp_tibble
          #    validObject(.Object)
              .Object
          }
)

# show method for TargetpResult object
setMethod("show",
          signature = 'TargetpResult',
          definition = function(object){
              cat(paste("An object of class", class(object), "\n"))
              if (length(object@seqList) != 0) {
                  print(elementNROWS(object@seqList))
              } else {
                  cat("all fasta slots are empty", '\n')          
              }
              cat('TargetP tabular output:', '\n')
              print(object@tp_tibble)
              
          })


#' Accessors for TargetpResult objects
#' @param theObject an object of TargetpResult class
#' @param tp_tibble parsed targetp output in tabular format
#' @export
#' @return TargetpResult object
#' @docType methods
#' @rdname TargetpResult_methods

setGeneric(
    name = "setTPtibble",
    def = function(theObject, tp_tibble)
    {
        standardGeneric("setTPtibble")
    }
)

#' @export
#' @rdname  TargetpResult_methods
#' @aliases setTPtibble

setMethod(
    f = "setTPtibble",
    signature = "TargetpResult",
    definition = function(theObject, tp_tibble)
    {
        theObject@tp_tibble <- tp_tibble
        validObject(theObject)
        return(theObject)
    }
)

setGeneric(
    name = "getTPtibble",
    def = function(theObject)
    {
        standardGeneric("getTPtibble")
    }
)

#' @export
#' @rdname  TargetpResult_methods
#' @aliases getTPtibble

setMethod(
    f = "getTPtibble",
    signature = "TargetpResult",
    definition = function(theObject)
    {
        return(theObject@tp_tibble)
    }
)
