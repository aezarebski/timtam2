library(dplyr)
library(reshape2)
library(ggplot2)
library(magrittr)
library(purrr)

out_dir <- "out"
bdsky_log <- file.path(out_dir, "ex2-bdsky-serial.log")
timtam_log <- file.path(out_dir, "ex2-timtam.log")
birth_rates <- c(4.5, 2.5) # 4.5 then 2.5 changing at time 2 (forwards)
birth_rate_change_time <- c(2.0)
death_rate <- 1.0
sampling_rate <- 0.5
occurrence_rate <- 0.5
final_prevalence <- 749
duration <- 5.0

# https://colorbrewer2.org/?type=diverging&scheme=PRGn&n=5
hex_codes <- c('#7b3294', # purple
               '#c2a5cf',
               '#f7f7f7',
               '#a6dba0',
               '#008837') # green


read_beast2_log <- function(filename) {
  read.csv(filename, sep = "\t", comment.char = "#")
}


bdsky_mcmc <- bdsky_log |>
  read_beast2_log() |>
  select(matches("birthRate*"))


timtam_mcmc <- timtam_log |>
  read_beast2_log() |>
  select(matches("birthRate*"), matches("*prevalence*")) |>
  mutate(
    prev_nb_mean = TimTam.prevalence.mean,
    prev_nb_size = ((TimTam.prevalence.mean ^ 2) /
                    (TimTam.prevalence.variance - TimTam.prevalence.mean))
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


## TODO Obviously everything below here needs to be fixed but this demonstrates
## the type of figure that should be created.

step_times <- c(0, birth_rate_change_time-1e-6, birth_rate_change_time + 1e-6, duration)

## this is a summary of the estimated birth rates from BDSky
bd_br_1 <- quantile(bdsky_mcmc$birthRate_BDSKY_Serial.1, probs=c(0.025, 0.5, 0.975))
bd_br_2 <- quantile(bdsky_mcmc$birthRate_BDSKY_Serial.2, probs=c(0.025, 0.5, 0.975))
tmp_bdsky <- data.frame(
  t = step_times,
  y = c(bd_br_1[2], bd_br_1[2], bd_br_2[2], bd_br_2[2]),
  ymin = c(bd_br_1[1], bd_br_1[1], bd_br_2[1], bd_br_2[1]),
  ymax = c(bd_br_1[3], bd_br_1[3], bd_br_2[3], bd_br_2[3])
)

tt_br_1 <- quantile(timtam_mcmc$birthRate.1, probs=c(0.025, 0.5, 0.975))
tt_br_2 <- quantile(timtam_mcmc$birthRate.2, probs=c(0.025, 0.5, 0.975))
tmp_timtam <- data.frame(
  t = step_times,
  y = c(tt_br_1[2], tt_br_1[2], tt_br_2[2], tt_br_2[2]),
  ymin = c(tt_br_1[1], tt_br_1[1], tt_br_2[1], tt_br_2[1]),
  ymax = c(tt_br_1[3], tt_br_1[3], tt_br_2[3], tt_br_2[3])
)

## This is a representation of the true birth rates used in the simulation.
tmp_sim <- data.frame(
  t = step_times,
  y = c(birth_rates[1], birth_rates[1], birth_rates[2], birth_rates[2])
)

tmp_gg <- ggplot() +
  geom_ribbon(
    data = tmp_bdsky,
    mapping = aes(x = t, ymin = ymin, ymax = ymax),
    alpha = 0.2,
    fill = hex_codes[4]
  ) +
  geom_line(
    data = tmp_bdsky,
    mapping = aes(x = t, y = y),
    colour = hex_codes[5]
  ) +
  geom_ribbon(
    data = tmp_timtam,
    mapping = aes(x = t, ymin = ymin, ymax = ymax),
    alpha = 0.2,
    fill = hex_codes[2]
  ) +
  geom_line(
    data = tmp_timtam,
    mapping = aes(x = t, y = y),
    colour = hex_codes[1]
    ) +
  geom_line(
    data = tmp_sim,
    mapping = aes(x = t, y = y),
    linetype = "dotted",
    size = 1
  ) +
  labs(x = "(Forwards) Time", y = "Birth Rate")

ggsave(filename = file.path(out_dir, "demonstration-plot.png"),
       plot = tmp_gg,
       height = 10.5, width = 14.8,
       ## A5 height = 14.8, width = 21.0,
       ## A6 height = 10.5, width = 14.8,
       ## A7 height = 7.4, width = 10.5,
       units = "cm")
