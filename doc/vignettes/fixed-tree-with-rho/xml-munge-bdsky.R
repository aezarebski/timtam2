library(xml2)

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

input_times_csv <- file.path(out_dir, "ape-sim-event-times.csv")
input_occ_csv <- file.path(out_dir, "occurrence-times.txt")
input_newick <- file.path(out_dir, "ape-sim-reconstructed-tree.newick")
simulation_duration <- duration
input_xml <- "bdsky-serial-2022-03-28.xml"
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

## ==============================================================================

#' Once we have the serially sampled analysis, we can make some further edits to
#' incorporate the contemporaneous sample into the specification.
rho_param <- list(
  id = "rho_BDSKY_Contemp.t:sequences",
  init_value = 0.3,
  beta_prior_alpha = 4.0,
  beta_prior_beta = 6.0
)

rho_param_node <- read_xml("<parameter spec='parameter.RealParameter' lower='0.0' name='stateNode' upper='1.0'></parameter>")
xml_set_attr(rho_param_node, attr = "id", value = rho_param$id)
xml_set_text(rho_param_node, value = as.character(rho_param$init_value))
xml_set_attr(
  xml_find_first(bdsky, "//distribution[@id='BDSKY_Serial.t:sequences']"),
  attr = "rho", value = paste0("@", rho_param$id)
)
xml_set_attr(
  xml_find_first(bdsky, "//distribution[@id='BDSKY_Serial.t:sequences']"),
  attr = "contemp", value = "true"
)

xml_add_child(xml_find_first(bdsky, "//state[@id='state']"), rho_param_node)
rho_prior_node <- read_xml(
  paste0(
    "<prior id='rhoPrior_BDSKY_Contemp.t:sequences' name='distribution'>",
    "<Beta id='Beta.0' name='distr'>",
    "<parameter id='RealParameter.1' spec='parameter.RealParameter' estimate='false' name='alpha'></parameter>",
    "<parameter id='RealParameter.2' spec='parameter.RealParameter' estimate='false' name='beta'></parameter>",
    "</Beta>",
    "</prior>"
    )
)
xml_set_attr(rho_prior_node, attr = "id", value = paste0("@", rho_param$id))
xml_set_text(
  xml_find_first(rho_prior_node, "//parameter[@name='alpha']"),
  value = as.character(rho_param$beta_prior_alpha)
)
xml_set_text(
  xml_find_first(rho_prior_node, "//parameter[@name='beta']"),
  value = as.character(rho_param$beta_prior_beta)
)

xml_add_child(
  xml_find_first(bdsky, "//distribution[@id='prior']"),
  rho_prior_node
)

rho_operator_node <- read_xml("<operator id='rhoScaler_BDSKY_Contemp.t:sequences' spec='ScaleOperator' weight='1.0' />")
xml_set_attr(
  rho_operator_node,
  attr = "parameter",
  value = paste0("@", rho_param$id)
)

xml_add_sibling(
  xml_find_first(bdsky, "//operator[@id='samplingProportionScaler_BDSKY_Serial.t:sequences']"),
  rho_operator_node
)

rho_log_node <- read_xml(sprintf("<log idref='%s' />", rho_param_id))
xml_add_child(xml_find_first(bdsky, "//logger[@id='tracelog']"), rho_log_node)

## ==============================================================================

xml_set_attr(
  xml_find_first(bdsky, "//logger[@id='tracelog']"),
  attr = "fileName",
  value = "bdsky-posterior.log")

xml_set_attr(
  xml_find_first(bdsky, "//run[@id='mcmc']"),
  attr = "chainLength",
  value = "100000"
)

## Write the result to file
write_xml(x = bdsky, file = output_again_xml)

## Format the result for easier reading.
system2(
  "js-beautify",
  args = c("-f", output_again_xml, "-r", "--type", "html"),
  wait = TRUE,
  timeout = 0)
