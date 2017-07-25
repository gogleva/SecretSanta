### SecretSanta workflow----
### Install signalPs (multiple versions if required) -----
### Put paths in the separate tab-delimeted file and upload them using manage_paths()----
secret_paths <- manage_paths("SecretSanta/inst/extdata/sample_paths") #this will leave in the environment
### Run signalp
result <- signalp(proteins = "SecretSanta/inst/extdata/sample_prot.fasta", organism_type = 'euk', version = 4)
result2 <- signalp(proteins = "SecretSanta/inst/extdata/sample_prot.fasta", organism_type = 'euk', version = 2)

### Parse signalP2/3 output:
res2 <- parse_signalp("SecretSanta/inst/extdata/sample_prot_signalp2_out", input_type = "path")

### run tmhmm:

#tmhmm()


