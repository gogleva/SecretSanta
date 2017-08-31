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
- TargetP
- SP cleavage
- TMHMM on mature peptides
- (K/H)DEL removal
- WoLFpsort
- M-slicer
- pipe them

2) annotate
- InterProSacan
- HMMs for specific motifs:
    - RxLR
    - WY
    - ER
    - CRN
    - NLS signals

3) Compare annotations and sequence-based prediction of signal peptides -> remove false positives

4) visualise inter-species comparisons


Expansion:
integrate some/all of the following tools:
http://bioinformatics.ysu.edu/tools/subcell.html


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
- -- signalp2
- -- signalp3
- -- signalp4 or 4.1
- -- targetp
- -- wolfpsort
- -- TMHMM
- clean/complete all function documentation
- unit tests
- targetp keeps failing; solution: shorten the path!
- exsaustive unit test of inputs/outputs
- function and class to run targetP
- function: M-slicer - for alternative translation start sites (could be mitochondrial?)
- speed up signalp run by splitting large input file into several smaller ones and run them as a series of small jobs
- Add check for empty inputs/ Message if the tool returns 0 candidates to pass further
- replace simple signalp with parallelised version
- unit tests and integration tests for parallelised signalp
- pesky closed connection warnings fixed
- clean tmp files after the TMHMM run (solved)
- does parLapply work on Windows? what to use if not? (they say it does), even if not -
stand alone CBS tools are available only for Unix, se we should not care about Windows so much.
- parallel targetp: replace simple targetp, microbenchmark, test on cluster, run unit and integration tests
- parallel TMHMM (implemented, but slower than simple tmhmm)
- parallel wolfpsort (not required, is already fast enough)
- consider having separate helper functions that are used by all/several runners (solved):
    - truncate_seq
    - smth else?
- add progress bars for parallelised tools, may be -> creates performance overhead - so, dropped    
- Write vignettes for SecretSanta predict
    - Installation instructions (done)
    - Add links and references to the tool papers in function docs and vignette (done)
- Vignettes: individual functions:
    - signalp (done)
    - targetp (done)
    - TMHMM (done)
    - WoLFpsort (done)
    - check_khdel (done)
    - M-slicer:slice (done)
    - M-slicer: rescue still fails with signalp (done)
- R CMD check: fix errors
- R CMD check: fix warnings    
- R CMD build: no warnings or errors
- R CMD check: fix notes
- Bio check: fix errors
- Bio check: fix warnings
- Vignettes: pipelines:
    - summary on starters and pipers;
    - short pipeline
- vignette: add types of organisms in the tool intros    
- external dependencies:is there a better way than manage_paths?
- signalp: weird Cmax scores -  check
- signalp: separate check for invalid organism
- signalp: separate check for invalid verison

TO DO (required minimum):

- parallel runners: Biocheck might complain about library exports
- add unit tests for large inputs in signalp - store them in a separate repo,
take too long to run for R CMD check:
    - 1K,
    - 2K,
    - long sequences
- why parallel version of TMHMM is so much slower ---> do more profiling
- export cluster environment: try to get rid of my_pa (done for signalp, requires testing)
- code profiling to speed-up wrappers
- fix some ugly bits in the code ---> profiling again. Aim: not to be much slower that the wrapped tool
- Individual (and cool!) plotting functions for the outputs
- Biocheck: fix notes (not absolutely required)
- Clean function imports (::) when they are necessary, when are not?
- do not need TMHMM_parallel if it is so slow?
- prediction graphs (?) when each candidate was rejected -->?

### Installation:

for private git_hub repo:
```
install_github("gogleva/SecretSanta", auth_token = 'd0b6ce8fcc6e31ecb98bab6ef273ccab6eb4227a')
```

### Coverage estimation with covr:

```
install_github("jimhester/covr")
library(covr)
my_cov <- package_coverage("SecretSanta")
```

Results, 09.08.2107:

```
SecretSanta Coverage: 77.84%
R/archived_functions.R: 0.00%
R/objects.R: 19.23%
R/manage_paths.R: 73.91%
R/check_khdel.R: 95.00%
R/parse_signalp.R: 97.37%
R/run_tmhmm.R: 97.37%
R/run_signalp.R: 98.39%
R/m-slicer.R: 100.00%
R/run_targetp.R: 100.00%
R/run_wolfpsort.R: 100.00%
```
