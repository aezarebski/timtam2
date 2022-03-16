library(ggplot2)
library(cowplot)
library(jsonlite)
library(XML)

posterior_samples <- read.csv(
  "timtam-posterior-multi.log",
  sep = "\t",
  comment.char = "#")
true_parameters <- xmlToList(xmlParse("multi-nu-params.xml"))





posterior_samples$deathRate <- as.numeric(true_parameters$configuration$parameters["deathRate"])
true_parameters$rNaught <- as.numeric(true_parameters$configuration$parameters["birthRate"]) / (as.numeric(true_parameters$configuration$parameters["deathRate"]) + as.numeric(true_parameters$configuration$parameters["samplingRate"]) + as.numeric(true_parameters$configuration$parameters["occurrenceRate"]))

true_parameters$prevalence <- read_json("out-multi/ape-sim-final-prevalence.json",
                             simplifyVector = TRUE)

posterior_samples$prevalence <- rnbinom(
  n = nrow(posterior_samples),
  size = exp(posterior_samples$TimTam.prevalence.lnR),
  prob = 1 - exp(posterior_samples$TimTam.prevalence.lnP))
posterior_samples$approxOccurrenceRate <- posterior_samples$nuProb / 0.5
posterior_samples$rNaught <- posterior_samples$birthRate / (posterior_samples$deathRate + posterior_samples$samplingRate + posterior_samples$approxOccurrenceRate)

g1 <- ggplot() +
  geom_hex(data = posterior_samples,
           mapping = aes(x = rNaught, y = prevalence),
           bins = 20) +
  geom_point(data = data.frame(rNaught=true_parameters$rNaught,prevalence=true_parameters$prevalence),
             mapping = aes(x = rNaught, y = prevalence),
             size = 4,
             colour = "red") +
  labs(x = "R-naught", y = "Prevalence") +
  scale_fill_continuous(low = "white", high = "black") +
  theme_classic() +
  theme(legend.position = "none")

g2 <- ggplot() +
  geom_histogram(data = posterior_samples,
                 mapping = aes(x = approxOccurrenceRate, y = ..density..),
                 bins = 20) +
  geom_vline(xintercept = as.numeric(true_parameters$configuration$parameters["occurrenceRate"]), colour = "red", linetype = "dashed") +
  xlim(c(0, 1)) +
  labs(x = "Unsequenced sampling probability", y = NULL) +
  theme_classic()

ggsave(
  filename = "out-multi/posterior-plot-multi.png",
  plot = plot_grid(g1, g2),
  height = 10.5,
  width = 21.0,
  units = "cm"
)
