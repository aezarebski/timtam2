library(ape)
library(phangorn)
x <- read.tree("out-multi/ape-sim-reconstructed-tree.newick")
seqs <- simSeq(x, l = 1, rate = 0)
write.phyDat(seqs, file = "out-multi/sequences.fasta", format = "fasta")

x <- read.csv("out-multi/ape-sim-aggregated-event-times.csv")
duration <- 4.0
root_length <- abs(x[x$event == "origin", "time"])
disaster_times <- x[x$event == "nu", "time"] + root_length
disaster_counts <- x[x$event == "nu", "size"]
bwd_disaster_times <- duration - disaster_times

writeLines(text = c(paste(sprintf("%.2f", bwd_disaster_times),collapse = " "),
                    paste(sprintf("%d", disaster_counts),collapse = " ")),
           con = "out-multi/disasters.txt")
