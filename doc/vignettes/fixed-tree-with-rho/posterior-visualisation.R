library(ggplot2)
library(cowplot)
library(jsonlite)

posterior_samples <- read.csv("out/fixed-tree-with-rho.log",
                              sep = "\t", comment.char = "#")
true_parameters <- as.data.frame(read_json("my-params.json"))
true_parameters$rNaught <- true_parameters$birthRate / (true_parameters$deathRate + true_parameters$samplingRate + true_parameters$occurrenceRate)
true_parameters$prevalence <- read_json("out/ape-sim-final-prevalence.json",
                             simplifyVector = TRUE)

posterior_samples$prevalence <- rnbinom(
  n = nrow(posterior_samples),
  size = exp(posterior_samples$TimTam.prevalence.lnR),
  prob = 1 - exp(posterior_samples$TimTam.prevalence.lnP))

g <- ggplot() +
  geom_hex(data = posterior_samples,
           mapping = aes(x = birthRate / (deathRate + samplingRate + occurrenceRate), y = prevalence),
           bins = 20) +
  geom_point(data = true_parameters,
             mapping = aes(x = rNaught, y = prevalence),
             size = 4,
             colour = "red") +
  labs(x = "R-naught", y = "Prevalence") +
  scale_fill_continuous(low = "white", high = "black") +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  filename = "./out/posterior-plot.png",
  plot = g,
  height = 10.5,
  width = 10.5,
  units = "cm"
)
