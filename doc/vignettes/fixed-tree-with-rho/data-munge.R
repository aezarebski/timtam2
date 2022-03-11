library(ape)
library(phangorn)
set.seed(1)
x <- read.tree("out/ape-sim-reconstructed-tree.newick")
seqs <- simSeq(x, l = 1, rate = 0)
write.phyDat(seqs, file = "out/sequences.fasta", format = "fasta")

x <- read.csv("out/ape-sim-event-times.csv")
duration <- 4.0
root_length <- abs(x[x$event == "origin", "time"])
occurrence_times <- sort(x[x$event == "occurrence", "time"] + root_length)
bwd_occurrence_times <- duration - occurrence_times
writeLines(text = paste(bwd_occurrence_times,collapse = " "),
           con = "out/occurrence-times.txt")
