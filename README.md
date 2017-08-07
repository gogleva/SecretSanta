# SecretSanta
Species: Oomycetes, Fungi, (Plants)

- predict (works!)
- annotate
- remove false positives
- classify
- compare
- visualise

-----------
1) predict
- use signalP(2-3-4) functionality
- TargetP (?)
- SP cleavage
- TMHMM on mature peptides
- (K/H)DEL removal
- WoLFpsort
- M-slicer


Things to consider:

test functions

How to implement prediction:
- system calls for signalp(versions) installed locally
- current focus: pipe_signalp - function to combine versions for multi-step filtering


Expansion:
integrate some/all of the following tools:
http://bioinformatics.ysu.edu/tools/subcell.html


## Got it:

runner functions will return instances of output classes objects
but: should run differently based on class of the input object?
or can we have the same attributes?
will it cause confusion?

### things to consider:

SOLVED:
- shall we store intermediate fasta files? - Yes
- if so -  how? alternatively - create tmp dir and output intermediates there (solved: as tmp files)
- to do: test on real proteomes: ~30 K proteins
- Question/Problem: how to feed XStringSet object to system call: solved with tmp files (solved)
- Question/Problem1: output fasta file for gene_ids only present in output tibble (solved)
- Problem2: need a separate function to construct mature sequences based on sp_tibble and out_fasta (solved)
- alternative solution for protein localisation: wolfpsort - run locally, parse the output, etc (done)
- perform test runs with the supplied paths, check that all the tools are executable in principle (done)
  -- signalp2
  -- signalp3
  -- signalp4 or 4.1
  -- targetp
  -- wolfpsort
  -- TMHMM

- clean/complete all function documentation
- unit tests
- targetp keeps failing; Solutin: shorten the path!



TO DO:

- clean tmp files after the signalp run?
- speed up signalp run by splitting large input file into several smaller ones and run them as a series of small jobs
- plotting functions for individual peptides (will require full output) signalp
- prediction graphs (?) when each candidate was rejected -->?
- function and class to run targetP
- function: M-slicer - for alternative translation start sites (could be mitochondrial?)
- excessive unit test of inputs/outputs for pipers




### Installation:

for private git_hub repo:
```
install_github("gogleva/SecretSanta", auth_token = 'd0b6ce8fcc6e31ecb98bab6ef273ccab6eb4227a')
```
