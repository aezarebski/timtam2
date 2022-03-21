library(xml2)

config <- read_xml("my-params.xml")
opts <- xml_find_first(config, "//options")
params <- xml_find_first(config, "//parameters")
output_dir <- xml_attr(opts, "outputDirectory")
input_times_csv <- sprintf("%s/ape-sim-event-times.csv", output_dir)
input_occ_txt <- sprintf("%s/occurrence-times.txt", output_dir)
input_newick <- sprintf("%s/ape-sim-reconstructed-tree.newick", output_dir)
## input_xml <- "et-with-rho-2022-03-18.xml"
input_xml <- "et-with-rho-2022-03-21.xml"
output_xml <- gsub(
  pattern = ".xml",
  replacement = "-edited.xml",
  x = input_xml)

beauti_output <- read_xml(input_xml)

## Include the extra data corresponding to catastrophes and occurrences.
catast_times_node <- xml_find_first(beauti_output, "//distribution//catastropheTimes[@id='catastropheSchedule.t:ape-sim-sequences']")
xml_set_text(catast_times_node, value = "0.0")
occ_times_node <- xml_find_first(beauti_output, "//distribution//points[@id='occurrenceTimes.t:ape-sim-sequences']")
xml_set_text(occ_times_node, readLines(input_occ_txt))

## Remove from the state node the parameters that we are not interested in estimating.
xml_remove(xml_find_first(beauti_output, "//state//parameter[@id='occurrenceRate.t:ape-sim-sequences']"))
xml_remove(xml_find_first(beauti_output, "//state//parameter[@id='rootLength.t:ape-sim-sequences']"))

xml_remove(xml_find_first(beauti_output, "//distribution[@id='MultiMonophyleticConstraint.ape-sim-sequences']"))
xml_remove(xml_find_first(beauti_output, "//distribution//disasterTimes"))
xml_remove(xml_find_first(beauti_output, "//distribution//disasterCounts"))
xml_remove(xml_find_first(beauti_output, "//prior[@id='occurrenceRatePrior.t:ape-sim-sequences']"))
xml_set_attr(
  xml_find_first(beauti_output, "//distribution[@id='TimTam.t:ape-sim-sequences']"),
  attr = "omega",
  value = xml_attr(params, "occurrenceRate"))

xml_remove(xml_find_first(beauti_output, "//operator[@id='occurrenceRateScaler.t:ape-sim-sequences']"))
xml_remove(xml_find_first(beauti_output, "//operator[@id='occurrenceRateWalk.t:ape-sim-sequences']"))


xml_remove(xml_find_first(beauti_output, "//log[@idref='occurrenceRate.t:ape-sim-sequences']"))
xml_set_attr(
  xml_find_first(beauti_output, "//logger[@id='tracelog']"),
  attr = "fileName",
  value = "timtam-posterior.log")
xml_set_attr(
  xml_find_first(beauti_output, "//logger[@id='treelog.t:ape-sim-sequences']"),
  attr = "fileName",
  value = "timtam-posterior.trees")

## Set the rho parameter to the correct value...
xml_set_text(
  xml_find_first(beauti_output, "//parameter[@name='p']"),
  xml_attr(params, "rho"))

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
