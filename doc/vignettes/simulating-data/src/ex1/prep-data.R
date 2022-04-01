library(ape)
library(phangorn)
library(xml2)
library(jsonlite)
suppressPackageStartupMessages(library(dplyr))
set.seed(1)

ape_sim_xml <- "src/ex1/sim-params.xml"

config <- read_xml(ape_sim_xml)

sim_duration <- as.numeric(
  xml_attr(
    xml_find_first(config, "//parameters"),
    "duration")
)

out_dir <- xml_attr(
  xml_find_first(config, "//options"),
  "outputDirectory"
)


newick_path <- file.path(out_dir, "ape-sim-reconstructed-tree.newick")
recon_tree <- read.tree(newick_path)

## Relabel the tips so that the final tip is at time zero and time goes
## backwards. The times initially encoded on the tips go forwards with the
## origin as zero.
rt_abs_tip_times <- sapply(
  X = strsplit(recon_tree$tip.label, split = "_"),
  FUN = function(x) as.numeric(x[2])
)
rt_time_offset <- sim_duration - max(rt_abs_tip_times)
rt_tip_bwd_times <- sim_duration - (rt_abs_tip_times + rt_time_offset)
rt_tip_prefix <- sapply(
  X = strsplit(recon_tree$tip.label, split = "_"),
  FUN = function(x) x[1]
)
recon_tree$tip.label <- sprintf("%s_%.6f", rt_tip_prefix, rt_tip_bwd_times)

fasta_path <- file.path(out_dir, "sequences.fasta")
seq_sim <- simSeq(recon_tree, l = 1, rate = 0)
write.phyDat(
  seq_sim,
  file = fasta_path,
  format = "fasta"
)

event_times_df <- read.csv(file.path(out_dir, "ape-sim-event-times.csv"))

## Print out a table of the results.
print(
  rename(
    as.data.frame(table(event_times_df$event)),
    event = Var1,
    frequency = Freq
  )
)

root_length <-
  as.numeric(abs(select(filter(event_times_df, event == "origin"), time)))

abs_occ_times <-
  sort(filter(event_times_df, event == "occurrence")$time + root_length)

bwd_occ_times <- sim_duration - (abs_occ_times + rt_time_offset)

writeLines(text = paste(bwd_occ_times, collapse = " "),
           con = file.path(out_dir, "occurrence-times.txt"))
jsonlite::write_json(
            x = list(
              final_prevalence = read_json(
                "out/ex1/ape-sim-final-prevalence.json",
                simplifyVector = TRUE
              ),
              newick_tree = write.tree(recon_tree),
              root_length = root_length,
              fasta = readLines(fasta_path),
              backward_occurrence_times = bwd_occ_times
            ),
            path = file.path(out_dir, "ex1.json"),
            auto_unbox = T,
            pretty = TRUE,
            digits = 16
          )
