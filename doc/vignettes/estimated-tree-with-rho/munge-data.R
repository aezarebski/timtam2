library(xml2)

config <- read_xml("my-params.xml")
params <- xml_find_first(config, "//parameters")
opts <- xml_find_first(config, "//options")
output_dir <- xml_attr(opts, "outputDirectory")
duration <- as.numeric(xml_attr(params, "duration"))
x <- read.csv(sprintf("%s/ape-sim-event-times.csv", output_dir))
root_length <- abs(x[x$event == "origin", "time"])
occurrence_times <- sort(x[x$event == "occurrence", "time"] + root_length)
bwd_occurrence_times <- duration - occurrence_times
writeLines(text = paste(bwd_occurrence_times,collapse = " "),
           con = sprintf("%s/occurrence-times.txt", output_dir))
