### SecretSanta workflow----
### Install signalPs (multiple versions if required) -----
### Put paths in the separate tab-delimeted file and upload them using manage_paths()----
secret_paths <- manage_paths("SecretSanta/inst/extdata/sample_paths") #this will leave in the environment

### Example pipe pipeline:

#initialise SignalpResult object
inp <- SignalpResult()

#read fasta file in AAStringSet object
aa <- readAAStringSet("SecretSanta/inst/extdata/sample_prot.fasta", use.names = TRUE)

#assign this object to the input_fasta slot of SignalpResult object
inp <- setInfasta(inp, aa)

#run signalp2 on the initial file:
step1_sp2 <- signalp(inp, version = 2, 'euk', run_mode = "starter")

#run signalp3 on the result object, will automatically pass out_fasta slot to signalp3:
step2_sp3 <- signalp(step1_sp2, version = 3, 'euk', run_mode = "piper")

#run signalp4 on the result object, will automatically pass out_fasta slot to signalp4:
step3_sp4 <- signalp(step2_sp3, version = 4, 'euk', run_mode = "piper")

#run TMHMM on the previous step

step4_TMHMM <- tmhmm(step3_sp4)

#check termial ER retention signals

step5_ER <- check_khdel(step4_TMHMM, 'piper')

# predict subcellular localisation (wolf psort)


