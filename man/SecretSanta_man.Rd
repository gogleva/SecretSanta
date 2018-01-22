\title{SecretSanta, main functionality}
\name{SecretSanta}
\alias{SecretSanta}
\description{
The SecretSanta package provides an R interface for the integrative prediction of extracellular proteins that are secreted via classical pathways. 
}
\details{
  Main functionality includes:
  \describe{
    \item{\code{\strong{manage_paths()}}}{run tests with external dependencies to ensure correct installation;}
    \item{\code{\strong{signalp()}}}{predict signal peptides with SignalP 2.0, SignalP 3.0 or SignalP 4.1;}
    \item{\code{\strong{tmhmm()}}}{predict transmembrane domains with TMHMM 2.0;}
    \item{\code{\strong{topcons()}}}{parse predictions of transmemrane domains performed by TOPCONS2;}
    \item{\code{\strong{targetp()}}}{predict subcellular localisation with TargetP 1.1;}
    \item{\code{\strong{wolfpsort()}}}{predict subcellular localisation with WoLF PSORT;}
    \item{\code{\strong{check_khdel()}}}{check C-terminal ER-retention signals;}
    \item{\code{\strong{m_slicer()}}}{generate proteins with alternative translation start sites;}
    \item{\code{\strong{ask_uniprot()}}}{fetch known subcellular location data from UniprotKB based on uniprot ids.}
  }
}
