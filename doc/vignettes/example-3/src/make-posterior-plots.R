library(dplyr)
library(reshape2)
library(ggplot2)
library(magrittr)
library(purrr)
library(jsonlite)

out_dir <- "out"
timtam_log <- file.path(out_dir, "ex3-timtam.log")
true_birth_rate <- 2.5
death_rate <- 1.0
sampling_rate <- 1.0
final_prevalence <- read_json(file.path("data", "ex3.json"))$final_prevalence
duration <- 10.0

# https://colorbrewer2.org/?type=diverging&scheme=PRGn&n=5
hex_codes <- c('#7b3294', # purple
               '#c2a5cf',
               '#f7f7f7',
               '#a6dba0',
               '#008837') # green


read_beast2_log <- function(filename) {
  read.csv(filename, sep = "\t", comment.char = "#")
}



timtam_mcmc <- timtam_log |>
  read_beast2_log() |>
  tail(250) |>
  select(matches("birthRate*"), matches("*prevalence*")) |>
  rename(
    birth_rate = birthRate,
    prevalence.mean = BDSKY_Serial.prevalence.mean,
    prevalence.variance = BDSKY_Serial.prevalence.variance
  ) |>
  mutate(
    prev_nb_mean = prevalence.mean,
    prev_nb_size = ((prevalence.mean ^ 2) /
                    (prevalence.variance - prevalence.mean))
  )

## Here we make a figure showing the estimate of the prevalence at the end of the simulation.
prev_df <- data.frame(
    prevalence = rnbinom(
      n = nrow(timtam_mcmc),
      mu = timtam_mcmc$prev_nb_mean,
      size = timtam_mcmc$prev_nb_size
    )
)

prev_gg <- ggplot(
  data = prev_df,
  mapping = aes(x = prevalence, y = ..density..)
) +
  geom_histogram(colour = hex_codes[1], fill = hex_codes[2], alpha=0.2) +
  geom_vline(xintercept = final_prevalence, linetype = "dotted", size = 1) +
  labs(x = "Prevalence", y = "Density")

ggsave(filename = file.path(out_dir, "prevalence-estimate.png"),
       plot = prev_gg,
       height = 10.5, width = 14.8,
       units = "cm")

## Here we make a figure showing the posterior distribution of the birth rate.
br_gg <- ggplot() +
  geom_histogram(
    data = timtam_mcmc,
    mapping = aes(x = birth_rate, y = ..density..),
    colour = hex_codes[1], fill = hex_codes[2], alpha = 0.2
  ) +
  geom_vline(
    xintercept = true_birth_rate,
    linetype = "dotted", size = 1) +
  stat_function(
    fun = function(x) dgamma(x = x, shape = 5, rate = 2.0),
    geom = "line",
    colour = hex_codes[1],
    size = 1
  ) +
  labs(x = "Birth rate", y = "Density")

ggsave(filename = file.path(out_dir, "birth-rate-estimate.png"),
       plot = br_gg,
       height = 10.5, width = 14.8,
       units = "cm")
