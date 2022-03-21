library(ggplot2)
library(cowplot)
library(jsonlite)
library(xml2)

posterior_samples <- read.csv("timtam-posterior.log",
                              sep = "\t", comment.char = "#")
## posterior_samples <- read.csv("out/estimated-tree-with-rho-via-beauti.log",
##                               sep = "\t", comment.char = "#")

sim_xml <- read_xml("my-params.xml")
xml_params <- xml_find_first(sim_xml, "//parameters")
true_parameters <- list(
  deathRate = as.numeric(xml_attr(xml_params, attr = "deathRate")),
  birthRate = as.numeric(xml_attr(xml_params, attr = "birthRate")),
  samplingRate = as.numeric(xml_attr(xml_params, attr = "samplingRate")),
  occurrenceRate = as.numeric(xml_attr(xml_params, attr = "occurrenceRate"))
)
posterior_samples$deathRate <- true_parameters$deathRate
true_parameters$rNaught <- with(
  true_parameters,
  birthRate / (deathRate + samplingRate + occurrenceRate)
)
true_parameters$prevalence <- read_json("out/ape-sim-final-prevalence.json",
                                        simplifyVector = TRUE)
true_final_prev <- read_json("out/ape-sim-final-prevalence.json",
                             simplifyVector = TRUE)



prev_post_samples <- rnbinom(
  n = 10 * nrow(posterior_samples),
  size = exp(posterior_samples$TimTam.prevalence.lnR),
  prob = 1 - exp(posterior_samples$TimTam.prevalence.lnP))

g1 <- ggplot() +
  geom_hex(data = posterior_samples,
           mapping = aes(x = birthRate, y = samplingRate),
           bins = 20) +
  geom_point(data = as.data.frame(true_parameters),
             mapping = aes(x = birthRate, y = samplingRate),
             size = 4,
             colour = "red") +
  labs(x = "Birth rate", y = "Sampling rate") +
  scale_fill_continuous(low = "white", high = "black") +
  theme_classic() +
  theme(legend.position = "none")

g2 <- ggplot(mapping = aes(x = x, y = ..density..)) +
  geom_histogram(data = data.frame(x = prev_post_samples),
                 bins = 20) +
  geom_vline(xintercept = true_final_prev,
             colour = "red",
             size = 2,
             linetype = "dashed") +
  labs(x = "Final prevalence", y = "Posterior distribution") +
  theme_classic()

g3 <- ggplot(mapping = aes(x = birthRate / (samplingRate + deathRate + true_parameters$occurrenceRate),
                           y = ..density..)) +
  geom_histogram(data = posterior_samples,
                 bins = 20) +
  geom_vline(xintercept = true_parameters$rNaught,
             colour = "red",
             size = 2,
             linetype = "dashed") +
  labs(x = "R-naught", y = "Posterior distribution") +
  theme_classic()

g <- plot_grid(g3, g2)


ggsave(
  filename = "./out/posterior-plot.png",
  plot = g,
  height = 10.5,
  width = 22.5,
  units = "cm"
)
