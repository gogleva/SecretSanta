### SecretSanta workflow----
### Install signalPs (multiple versions if required) -----
### Put paths in the separate tab-delimeted file and upload them using manage_paths()----
secret_paths <- manage_paths("SecretSanta/inst/extdata/sample_paths") #this will leave in the environment
### Run signalp
result <- signalp(proteins = "SecretSanta/inst/extdata/sample_prot.fasta", organism_type = 'euk', version = 4)






