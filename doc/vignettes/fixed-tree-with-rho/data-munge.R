library(ape)
library(phangorn)
library(xml2)
set.seed(1)

config <- read_xml("sim-const-params.xml")
duration <- as.numeric(
  xml_attr(
    xml_find_first(config, "//parameters"),
    "duration")
)
out_dir <- xml_attr(
  xml_find_first(config, "//options"),
  "outputDirectory"
)

tree <- read.tree(file.path(out_dir, "ape-sim-reconstructed-tree.newick"))
write.phyDat(
  simSeq(tree, l = 1, rate = 0),
  file = file.path(out_dir, "sequences.fasta"),
  format = "fasta"
)

x <- read.csv(file.path(out_dir, "ape-sim-event-times.csv"))
root_length <- abs(x[x$event == "origin", "time"])
occurrence_times <- sort(x[x$event == "occurrence", "time"] + root_length)
bwd_occurrence_times <- duration - occurrence_times
writeLines(text = paste(bwd_occurrence_times, collapse = " "),
           con = file.path(out_dir, "occurrence-times.txt"))
