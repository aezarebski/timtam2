library(ape)
library(phangorn)
library(xml2)
library(jsonlite)
suppressPackageStartupMessages(library(dplyr))
set.seed(1)

ape_sim_xml <- "src/ex1/sim-params.xml"

config <- read_xml(ape_sim_xml)

duration <- as.numeric(
  xml_attr(
    xml_find_first(config, "//parameters"),
    "duration")
)

out_dir <- xml_attr(
  xml_find_first(config, "//options"),
  "outputDirectory"
)

newick_path <- file.path(out_dir, "ape-sim-reconstructed-tree.newick")
tree <- read.tree(newick_path)

fasta_path <- file.path(out_dir, "sequences.fasta")
write.phyDat(
  simSeq(tree, l = 1, rate = 0),
  file = fasta_path,
  format = "fasta"
)

x <- read.csv(file.path(out_dir, "ape-sim-event-times.csv"))

## Print out a table of the results.
print(
  rename(
    as.data.frame(table(x$event)),
    event = Var1,
    frequency = Freq
  )
)

root_length <- as.numeric(abs(select(filter(x, event == "origin"), time)))
occurrence_times <- sort(x[x$event == "occurrence", "time"] + root_length)
bwd_occurrence_times <- duration - occurrence_times
writeLines(text = paste(bwd_occurrence_times, collapse = " "),
           con = file.path(out_dir, "occurrence-times.txt"))
jsonlite::write_json(
            x = list(
              final_prevalence = read_json(
                "out/ex1/ape-sim-final-prevalence.json",
                simplifyVector = TRUE
              ),
              newick_tree = readLines(newick_path),
              root_length = root_length,
              fasta = readLines(fasta_path),
              backward_occurrence_times = bwd_occurrence_times
            ),
            path = file.path(out_dir, "ex1.json"),
            auto_unbox = T,
            pretty = TRUE,
            digits = 16
          )
