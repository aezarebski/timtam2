library(phangorn)

alignment_fasta <- "out/ape-sim-sequences.fasta"
tree_nexus <- "out/phangorn-nj-tree.nexus"

alignment <- read.phyDat(
  file = alignment_fasta,
  format = "fasta",
  type = "DNA"
)
unrooted_tree <- nj(as.DNAbin(alignment))

write.nexus(unrooted_tree, file = tree_nexus)
