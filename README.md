## 1. Background
**SecretSanta** provides an R interface aiding integrative prediction of
extracellular proteins secreted via classical secretion pathway, i.e. requiring
presence of a signal peptide.

Secretome prediction often involves multiple steps. Typically it starts with
prediction of short signal peptides at the N-terminus end of a protein
(**Figure 1, a**). Next, it is crucial to ensure the absence of motifs and
domains preventing the protein from being secreted despite the presence of the
signal peptide. These sequences include transmembrane domains, short ER lumen
retention signals, mitochondria and plastid targeting signals
(**Figure 1, b-d**). The ultimate aim of a secretome prediction pipeline is to
pick secreted proteins shown in **Figure 1, a** and filer out
those shown in **Figure 1, b-d**.


\n

A number of excellent command line tools and web-interfaces exist to perform
prediction of individual motifs and domains
([signalp](http://www.cbs.dtu.dk/services/SignalP/),
[targetp](http://www.cbs.dtu.dk/services/TargetP/),
[TMHMM](http://www.cbs.dtu.dk/services/TMHMM/),
[WolfPsort](https://github.com/fmaguire/WoLFPSort)), however the interface
allowing to combine the outputs in a single flexible workflow is lacking.

\n

**SecretSanta** package attempts to bridge this gap. It provides wrapper
functions around existing command line tools for prediction of signal peptides
and protein subcellular localisation. The wrappers are designed to work
together by producing standardized output as an instance of ``CBSResult`` superclass.
