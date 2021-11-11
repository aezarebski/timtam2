library(ggplot2)

posterior_samples <- read.csv("test.1636648929209.log", sep = "\t")
true_parameters <- as.data.frame(jsonlite::read_json("my-params.json"))

g <- ggplot() +
  geom_hex(data = posterior_samples,
           mapping = aes(x = myLambda, y = myPsi),
           bins = 20) +
  geom_point(data = true_parameters,
             mapping = aes(x = birthRate, y = samplingRate),
             size = 4,
             colour = "red") +
  labs(x = "Birth rate", y = "Sampling rate") +
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
