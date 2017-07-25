# SecretSanta
Species: Oomycetes, Fungi, (Plants)

- predict
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


Things to consider:

test functions

How to implement prediction:
- system calls for signalp(versions) installed locally
- current focus: pipe_signalp - function to combine versions for multi-step filtering


Expansion:
integrate some/all of the following tools:
http://bioinformatics.ysu.edu/tools/subcell.html


## Got it:

runner fucntions will return instances of output classes objects
but: should run differently based on class of the input object? 
or can we have the same attributes?
will it cause confusion?

### things to consider:

-- shall we store intermediate fast files?
-- if so -  how? alternatively - create tmp dir and output intermediates there
-- to do: test on real proteomes: ~30 K proteins


