#!/usr/bin/env Rscript

library(RSVSim)

path_output = "/proj/breedmap/NOBACKUP/sim_data"
path_genome = "/proj/breedmap/REF/ARS-UCD1.2_Btau5.0.1Y.fa.gz"

delsizes = sample(100:10000,20,replace = TRUE) 
inssizes = sample(50:10000,40,replace = TRUE)
invsizes = sample(100:10000,27,replace = TRUE)

sim = simulateSV(output=path_output,genome = path_genome, dels = 20, sizeDels = delsizes,
		 ins = 40, sizeIns = inssizes, invs = 27, sizeInvs = invsizes, repeatBias = FALSE,
		 bpSeqSize = 39, seed = 777, verbose = FALSE)

metadata(sim)
 
