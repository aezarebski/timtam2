library(xml2)

input_times_csv <- "out/ape-sim-event-times.csv"
input_occ_csv <- "out/occurrence-times.txt"
input_newick <- "out/ape-sim-reconstructed-tree.newick"
simulation_duration <- "4.0"
input_xml <- "bdsky-serial-2022-03-14.xml"
output_xml <- gsub(
  pattern = ".xml",
  replacement = "-edited.xml",
  x = input_xml)
output_again_xml <- gsub(
  pattern = "-serial",
  replacement = "",
  x = output_xml)

bdsky <- read_xml(input_xml)

xml_remove(xml_find_first(bdsky, "//tree"))
xml_remove(xml_find_first(bdsky, "//parameter[@id='clockRate.c:sequences']"))
xml_remove(xml_find_first(bdsky, "//parameter[@id='origin_BDSKY_Serial.t:sequences']"))
xml_set_attr(xml_find_first(bdsky, "//parameter[@id='reproductiveNumber_BDSKY_Serial.t:sequences']"), attr = "dimension", value = NULL)
xml_set_text(xml_find_first(bdsky, "//parameter[@id='becomeUninfectiousRate_BDSKY_Serial.t:sequences']"), "3.0")

xml_remove(xml_find_first(bdsky, "//distribution[@id='MultiMonophyleticConstraint.sequences']"))

## Set the starting tree to be the reconstructed tree.
fixed_tree_node <- read_xml("<init spec='beast.util.TreeParser' taxa='@sequences' id='Tree.t:sequences' IsLabelledNewick='true' adjustTipHeights='false' />")
xml_set_attr(
  x = fixed_tree_node,
  attr = "newick",
  value = readLines(input_newick))
rand_tree_node <- xml_find_first(bdsky, "//init[@id='RandomTree.t:sequences']")
xml_replace(rand_tree_node, fixed_tree_node)


xml_set_attr(xml_find_first(bdsky, "//distribution[@id='BDSKY_Serial.t:sequences']"), attr = "origin", value = simulation_duration)

## Remove unused operators
xml_remove(xml_find_first(bdsky, "//operator[@id='BDSKY_SerialTreeScaler.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='BDSKY_SerialTreeRootScaler.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='BDSKY_SerialUniformOperator.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='BDSKY_SerialSubtreeSlide.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='BDSKY_SerialNarrow.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='BDSKY_SerialWide.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='BDSKY_SerialWilsonBalding.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='StrictClockRateScaler.c:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='strictClockUpDownOperator.c:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='deathRateScaler.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='deathRateWalk.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='rootLengthScaler.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='rootLengthWalk.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='nuProbWalk.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//operator[@id='origScaler_BDSKY_Serial.t:sequences']"))

xml_remove(xml_find_first(bdsky, "//log[@idref='likelihood']"))
xml_remove(xml_find_first(bdsky, "//log[@idref='treeLikelihood.sequences']"))
xml_remove(xml_find_first(bdsky, "//log[@id='TreeHeight.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//log[@idref='clockRate.c:sequences']"))
xml_remove(xml_find_first(bdsky, "//log[@idref='likelihood']"))
xml_remove(xml_find_first(bdsky, "//log[@idref='origin_BDSKY_Serial.t:sequences']"))
xml_remove(xml_find_first(bdsky, "//logger[@id='treelog.t:sequences']"))

xml_remove(xml_find_first(bdsky, "//distribution[@id='likelihood']"))

xml_remove(xml_find_first(bdsky, "//prior[@id='ClockPrior.c:sequences']"))
xml_remove(xml_find_first(bdsky, "//prior[@id='originPrior_BDSKY_Serial.t:sequences']"))




## Write the result to file
write_xml(x = bdsky, file = output_xml)

## Format the result for easier reading.
system2(
  "js-beautify",
  args = c("-f", output_xml, "-r", "--type", "html"),
  wait = TRUE,
  timeout = 0)


#' Once we have the serially sampled analysis, we can make some further edits to
#' incorporate the contemporaneous sample.

rho_param_id <- "rho_BDSKY_Contemp.t:sequences"
rho_param_node <- read_xml(sprintf("<parameter id='%s' spec='parameter.RealParameter' lower='0.0' name='stateNode' upper='1.0'>0.3</parameter>", rho_param_id))
xml_add_child(xml_find_first(bdsky, "//state[@id='state']"), rho_param_node)

xml_set_attr(
  xml_find_first(bdsky, "//distribution[@id='BDSKY_Serial.t:sequences']"),
  attr = "contemp",
  value = "true"
)
xml_set_attr(
  xml_find_first(bdsky, "//distribution[@id='BDSKY_Serial.t:sequences']"),
  attr = "rho",
  value = sprintf("@%s", rho_param_id)
)

rho_prior_node <- read_xml(paste0(c(sprintf("<prior id='rhoPrior_BDSKY_Contemp.t:sequences' name='distribution' x='@%s'>", rho_param_id),
  "<Beta id='Beta.0' name='distr'>",
  "<parameter id='RealParameter.1' spec='parameter.RealParameter' estimate='false' name='alpha'>4.0</parameter>",
  "<parameter id='RealParameter.2' spec='parameter.RealParameter' estimate='false' name='beta'>4.0</parameter>",
  "</Beta>",
  "</prior>"), collapse = ""))

xml_add_child(xml_find_first(bdsky, "//distribution[@id='prior']"), rho_prior_node)

rho_operator_node <- read_xml(sprintf("<operator id='rhoScaler_BDSKY_Contemp.t:sequences' spec='ScaleOperator' parameter='@%s' weight='1.0' />", rho_param_id))
xml_add_sibling(xml_find_first(bdsky, "//operator[@id='samplingProportionScaler_BDSKY_Serial.t:sequences']"), rho_operator_node)

rho_log_node <- read_xml(sprintf("<log idref='%s' />", rho_param_id))
xml_add_child(xml_find_first(bdsky, "//logger[@id='tracelog']"), rho_log_node)

xml_set_attr(
  xml_find_first(bdsky, "//logger[@id='tracelog']"),
  attr = "fileName",
  value = "bdsky-posterior.log")

xml_set_attr(
  xml_find_first(bdsky, "//run[@id='mcmc']"),
  attr = "chainLength",
  value = "1000000"
)

## Write the result to file
write_xml(x = bdsky, file = output_again_xml)

## Format the result for easier reading.
system2(
  "js-beautify",
  args = c("-f", output_again_xml, "-r", "--type", "html"),
  wait = TRUE,
  timeout = 0)
