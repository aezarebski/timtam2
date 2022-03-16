library(xml2)

input_times_csv <- "out/ape-sim-event-times.csv"
input_occ_csv <- "out/occurrence-times.txt"
input_newick <- "out/ape-sim-reconstructed-tree.newick"
input_xml <- "ft-with-nu-2022-03-15.xml"
output_xml <- gsub(
  pattern = ".xml",
  replacement = "-edited.xml",
  x = input_xml)

beauti_output <- read_xml(input_xml)

## Remove some nodes that we don't want.

unwanted_state_nodes <-paste(
  c("//state//tree",
    sprintf(
      "//state//parameter[@id='%s']",
      c(
        "clockRate.c:sequences",
        "rootLength.t:sequences",
        "Tree.t:sequences"))),
  collapse = "|")

xml_remove(
  xml_find_all(
    x = beauti_output,
    xpath = unwanted_state_nodes),
  free = TRUE)

## Adjust the initial value of the birth rate.
state_birth_rate <- xml_find_first(beauti_output, "//state//parameter[@id='birthRate.t:sequences']")
xml_text(state_birth_rate) <- "3.0"

## Adjust the occurrence times data.
occurrence_times <- xml_find_first(beauti_output, "//distribution//points[@id='occurrenceTimes.t:sequences']")
xml_text(occurrence_times) <- readLines(input_occ_csv)

## Adjust the disaster times schedule relative to the last sequenced sample.
tmp <- read.csv(input_times_csv)
last_seq_time <- max(tmp[tmp$event == "sampling","time"])
disaster_size <- nrow(tmp[tmp$event == "nu",])
disaster_time <- unique(tmp[tmp$event == "nu","time"])
disaster_schedule <- xml_find_first(
  beauti_output,
  "//distribution//disasterTimes[@id='disasterSchedule.t:sequences']")
xml_text(disaster_schedule) <- sprintf("%.7f", last_seq_time - disaster_time)
disaster_counts <- xml_find_first(
  beauti_output,
  "//distribution//disasterCounts[@id='disasterCounts.t:sequences']")
xml_text(disaster_counts) <- as.character(disaster_size)
rm(tmp)

## Set the starting tree to be the reconstructed tree.
fixed_tree_node <- read_xml("<init spec='beast.util.TreeParser' taxa='@sequences' id='Tree.t:sequences' IsLabelledNewick='true' adjustTipHeights='false' />")
xml_set_attr(
  x = fixed_tree_node,
  attr = "newick",
  value = readLines(input_newick))
rand_tree_node <- xml_find_first(beauti_output, "//init[@id='RandomTree.t:sequences']")
xml_replace(rand_tree_node, fixed_tree_node)

## Insert the root length properly and remove unused variables.
tmp <- read.csv(input_times_csv, nrows=1)
root_length <- abs(tmp[tmp$event == "origin","time"])
rm(tmp)
timtam_dist_node <- xml_find_first(beauti_output, "//distribution[@id='TimTam.t:sequences']")
xml_set_attr(
  timtam_dist_node,
  attr = "rootLength",
  value = root_length
)
## xml_set_attr(
##   timtam_dist_node,
##   attr = "nu",
##   value = NULL
## )
xml_set_attr(
  timtam_dist_node,
  attr = "mu",
  value = NULL
)

## Remove the some unused nodes from the distribution
xml_remove(xml_find_first(beauti_output, "//distribution[@id='MultiMonophyleticConstraint.sequences']"))
xml_remove(xml_find_first(beauti_output, "//distribution[@id='likelihood']"))

xml_remove(xml_find_first(beauti_output, "//prior[@id='ClockPrior.c:sequences']"))
xml_remove(xml_find_first(beauti_output, "//prior[@id='rootLengthPrior.t:sequences']"))
xml_remove(xml_find_first(beauti_output, "//operator[@id='TimTamTreeScaler.t:sequences']"))
xml_remove(xml_find_first(beauti_output, "//operator[@id='TimTamTreeRootScaler.t:sequences']"))
xml_remove(xml_find_first(beauti_output, "//operator[@id='TimTamUniformOperator.t:sequences']"))
xml_remove(xml_find_first(beauti_output, "//operator[@id='TimTamSubtreeSlide.t:sequences']"))
xml_remove(xml_find_first(beauti_output, "//operator[@id='TimTamNarrow.t:sequences']"))
xml_remove(xml_find_first(beauti_output, "//operator[@id='TimTamWide.t:sequences']"))
xml_remove(xml_find_first(beauti_output, "//operator[@id='TimTamWilsonBalding.t:sequences']"))
xml_remove(xml_find_first(beauti_output, "//operator[@id='StrictClockRateScaler.c:sequences']"))
xml_remove(xml_find_first(beauti_output, "//operator[@id='strictClockUpDownOperator.c:sequences']"))
xml_remove(xml_find_first(beauti_output, "//operator[@id='rootLengthScaler.t:sequences']"))
xml_remove(xml_find_first(beauti_output, "//operator[@id='rootLengthWalk.t:sequences']"))
xml_remove(xml_find_first(beauti_output, "//log[@idref='likelihood']"))
xml_remove(xml_find_first(beauti_output, "//log[@idref='treeLikelihood.sequences']"))
xml_remove(xml_find_first(beauti_output, "//log[@id='TreeHeight.t:sequences']"))
xml_remove(xml_find_first(beauti_output, "//log[@idref='clockRate.c:sequences']"))
xml_remove(xml_find_first(beauti_output, "//log[@idref='likelihood']"))
xml_remove(xml_find_first(beauti_output, "//log[@idref='rootLength.t:sequences']"))
xml_remove(xml_find_first(beauti_output, "//logger[@id='treelog.t:sequences']"))


xml_remove(xml_find_first(beauti_output, "//parameter[@name='p']"))
xml_remove(xml_find_first(beauti_output, "//catastropheTimes"))


xml_set_attr(
  xml_find_first(beauti_output, "//logger[@id='tracelog']"),
  attr = "fileName",
  value = "timtam-posterior.log")

xml_set_attr(
  xml_find_first(beauti_output, "//run[@id='mcmc']"),
  attr = "chainLength",
  value = "1000000"
)

## Write the result to file
write_xml(x = beauti_output, file = output_xml)

## Format the result for easier reading.
system2(
  "js-beautify",
  args = c("-f", output_xml, "-r", "--type", "html"),
  wait = TRUE,
  timeout = 0)
