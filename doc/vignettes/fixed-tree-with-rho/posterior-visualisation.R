library(ggplot2)
library(cowplot)
library(jsonlite)

posterior_samples <- read.csv("out/fixed-tree-with-rho.log",
                              sep = "\t", comment.char = "#")
true_parameters <- as.data.frame(read_json("my-params.json"))
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
  geom_point(data = true_parameters,
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

g <- plot_grid(g1, g2)


ggsave(
  filename = "./out/posterior-plot.png",
  plot = g,
  height = 10.5,
  width = 22.5,
  units = "cm"
)
