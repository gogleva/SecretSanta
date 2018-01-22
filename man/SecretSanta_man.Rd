\title{SecretSanta, main functionality}
\name{SecretSanta}
\alias{SecretSanta}
\description{
The SecretSanta package provides an R interface for the integrative prediction of extracellular proteins that are secreted via classical pathways. 
}
\details{
Main functionality includes:
\describe{
\item{\strong{\code{signalp()}}}{predict signal peptides with SignalP 2.0, SignalP 3.0 or SignalP 4.1;}
\item{\strong{\code{tmhmm()}}}{predict transmembrane domains with TMHMM 2.0;}
\item{\strong{\code{topcons()}}}{parse predictions of transmemrane domains performed by TOPCONS2;}
\item{\strong{\code{targetp()}}}{predict subcellular localisation with TargetP 1.1;}
\item{\strong{\code{wolfpsort()}}}{predict subcellular localisation with WoLF PSORT;}
\item{\strong{\code{check_khdel()}}}{check C-terminal ER-retention signals;}
\item{\strong{\code{m_slicer()}}}{generate proteins with alternative translation start sites;}
\item{\strong{\code{ask_uniprot()}}}{fetch known subcellular location data from UniprotKB based on uniprot ids.}
}
\cr
S4 classes and accessor methods:
\describe{
\item{\strong{CBSResult}}{Basic parent S4 class with \code{in_fasta} and \code{out_fasta} slots. Accessor functions:
\itemize{
\item{\code{getInfasta():} extract \code{in_fasta} slot;}
\item{\code{getOutfasta():} extract \code{out_fasta} slot;}
\item{\code{getFastas():} extract both/all fasta slots;}
}}
\item{\strong{SignalpResult}}{S4 class, organisig outputs of \code{signalp()} function. Inherits from \strong{CBSResult} class. Slots include: \code{in_fasta}, \code{out_fasta}, \code{mature_fasta}, \code{sp_tibble}, \code{sp_version}.    Specific accessor functions:
\itemize{
\item{\code{getMatfasta():} extract \code{mature_fasta} slot;}
\item{\code{getSPtibble():} extract SignalP summary stats;}
\item{\code{getSPversion():} extract version of SignalP used;}
}}
\item{\strong{TargetpResult}}{S4 class, organisig outputs of \code{targetp()} function. Inherits from \strong{CBSResult} class. Slots include: \code{in_fasta}, \code{out_fasta}, \code{tp_tibble}. Specific accessor function:
\itemize{
\item{\code{getTPtibble():} extract parsed TargetP output in tabular format;}
}}
\item{\strong{TMhmmResult}}{S4 class, organisig outputs of \code{tmhmm()} function. Inherits from \strong{CBSResult} class. Slots include: \code{in_fasta}, \code{out_fasta}, \code{in_mature_fasta}, \code{out_mature_fasta}, \code{tm_tibble}.    Specific accessor functions:
\itemize{
\item{\code{getInMatfasta():} extract \code{in_mature_fasta} slot;}
\item{\code{getOutMatfasta():} extract \code{out_mature_fasta} slot;}
\item{\code{getTMtibble():} extract TMHMM summary stats;}
}}
\item{\strong{TopconsResult}}{S4 class, organisig outputs of \code{topcons()} function. Inherits from \strong{CBSResult} class. Slots include: \code{in_fasta}, \code{out_fasta}, \code{top_tibble}. Specific accessor function:
\itemize{
\item{\code{getTOPtibble():} extract TOPCONS summary stats;}
}}
\item{\strong{WolfResult}}{S4 class, organisig outputs of \code{wolfpsort()} function. Inherits from \strong{CBSResult} class. Slots include: \code{in_fasta}, \code{out_fasta}, \code{wolf_tibble}. Specific accessor function:
\itemize{
\item{\code{getWOLFtibble():} extract WoLF PSORT summary;}
}}
\item{\strong{ErResult}}{S4 class, organisig outputs of \code{check_khdel()} function. Inherits from \strong{CBSResult} class. Slots include: \code{in_fasta}, \code{out_fasta}, \code{retained_fasta}. Specific accessor function:
\itemize{
\item{\code{getRetainedfasta():} extract \code{retained_fasta} slot.}
}}
}
}