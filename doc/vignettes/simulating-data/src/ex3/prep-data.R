library(ape)
library(xml2)
library(jsonlite)
suppressPackageStartupMessages(library(dplyr))
set.seed(1)

ape_sim_xml <- "src/ex3/sim-params.xml"

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

fasta_path <- file.path(out_dir, "ape-sim-sequences.fasta")
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


## We can make a NJ tree to use as a rough guess for what the resulting tree
## should look like and as a way to test if there is any signal in the data.
fasta_data <- read.FASTA(fasta_path)
seq_labels <- names(fasta_data)
seqs <- as.character(fasta_data)
num_seqs <- length(seqs)
dm <- matrix(0, num_seqs, num_seqs)
rownames(dm) <- seq_labels
colnames(dm) <- seq_labels
for (ix in 1:(num_seqs - 1)) {
  for (jx in ix:num_seqs) {
    jac <- sum(seqs[[ix]] != seqs[[jx]])
    dm[ix,jx] <- jac
    dm[jx,ix] <- jac
  }
}
nj_tree <- nj(dm)

## Finally we write all of this out to a blob.
writeLines(text = paste(bwd_occ_times, collapse = " "),
           con = file.path(out_dir, "occurrence-times.txt"))
jsonlite::write_json(
            x = list(
              final_prevalence = read_json(
                file.path(out_dir, "ape-sim-final-prevalence.json"),
                simplifyVector = TRUE
              ),
              fasta = readLines(fasta_path),
              nj = write.tree(phy = nj_tree),
              backward_occurrence_times = bwd_occ_times
            ),
            path = file.path(out_dir, "ex3.json"),
            auto_unbox = T,
            pretty = TRUE,
            digits = 16
          )
