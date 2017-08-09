### experimets with parallelisation

library(parallel)

# function to split tmp fasta into multiple fastas
# run signalp for each of them
# combine the result
# lapply is my friend

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

# large real-life file:
large_inp <- readAAStringSet(system.file("extdata", "Ppalm_prot_ALI_PLTG.fasta", package = "SecretSanta"))


parLapply()

stopCluster(cl)
