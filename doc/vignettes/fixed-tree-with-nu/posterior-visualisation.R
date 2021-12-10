library(ggplot2)
library(cowplot)
library(jsonlite)

posterior_samples <- read.csv("out/ft-with-nu-02.log",
                              sep = "\t", comment.char = "#")
true_parameters <- as.data.frame(read_json("my-params.json"))
posterior_samples$deathRate <- true_parameters$deathRate
true_parameters$rNaught <- true_parameters$birthRate / (true_parameters$deathRate + true_parameters$samplingRate + true_parameters$occurrenceRate)
true_parameters$prevalence <- read_json("out/ape-sim-final-prevalence.json",
                             simplifyVector = TRUE)

posterior_samples$prevalence <- rnbinom(
  n = nrow(posterior_samples),
  size = exp(posterior_samples$TimTam.prevalence.lnR),
  prob = 1 - exp(posterior_samples$TimTam.prevalence.lnP))
posterior_samples$rNaught <- posterior_samples$birthRate / (posterior_samples$deathRate + posterior_samples$samplingRate + posterior_samples$occurrenceRate)

g1 <- ggplot() +
  geom_hex(data = posterior_samples,
           mapping = aes(x = rNaught, y = prevalence),
           bins = 20) +
  geom_point(data = true_parameters,
             mapping = aes(x = rNaught, y = prevalence),
             size = 4,
             colour = "red") +
  labs(x = "R-naught", y = "Prevalence") +
  scale_fill_continuous(low = "white", high = "black") +
  theme_classic() +
  theme(legend.position = "none")

g2 <- ggplot() +
  geom_histogram(data = posterior_samples,
                 mapping = aes(x = nuProb, y = ..density..),
                 bins = 20) +
  stat_function(
    fun = function(x) dunif(x),
    geom = "line"
  ) +
  geom_vline(xintercept = 0.5, colour = "red", linetype = "dashed") +
  xlim(c(0, 1)) +
  labs(x = "Unsequenced sampling probability", y = NULL) +
  theme_classic()

ggsave(
  filename = "./out/posterior-plot.png",
  plot = plot_grid(g1, g2),
  height = 10.5,
  width = 21.0,
  units = "cm"
)
