library(ggplot2)
library(cowplot)
library(jsonlite)

## posterior_samples <- read.csv("out/fixed-tree-with-rho-01.log",
posterior_samples <- read.csv("out/timtam-posterior.log",
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

ggsave(
  filename = "./out/posterior-plot.png",
  plot = g1,
  height = 10.5,
  width = 10.5,
  units = "cm"
)

## bdsky_posterior_samples <- read.csv("out/fixed-tree-with-rho-03.log",
bdsky_posterior_samples <- read.csv("out/bdsky-posterior.log",
                                    sep = "\t", comment.char = "#")
plot_df <- data.frame(model = rep(c("timtam", "bdsky"), each = nrow(posterior_samples)),
                      rNaught = c(posterior_samples$rNaught, bdsky_posterior_samples$reproductiveNumber_BDSKY_Serial))

g2 <- ggplot() +
  geom_histogram(data = plot_df, mapping = aes(x = rNaught)) +
  geom_vline(xintercept = true_parameters$rNaught, linetype = "dashed") +
  labs(x = "R-naught", y = "Posterior samples") +
  facet_wrap(~model) +
  theme_classic()

ggsave(filename = "./out/r-naught-comparison.png",
       plot = g2,
       height = 14.8, width = 21.0,
       units = "cm")
