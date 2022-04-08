library(dplyr)
library(reshape2)
library(ggplot2)
library(magrittr)
library(purrr)
library(jsonlite)
library(stringr)
library(cowplot)
library(whisker)

## TODO these variables should be read out of the input data rather than
## hard-coded here.

out_dir <- "out"
bdsky_log <- file.path(out_dir, "ex1-bdsky-serial.log")
timtam_log <- file.path(out_dir, "ex1-timtam.log")
birth_rate <- 3.0
death_rate <- 1.0
sampling_rate <- 0.5
occurrence_rate <- 0.5
final_prevalence <- 1025


read_beast2_log <- function(filename) {
  read.csv(filename, sep = "\t", comment.char = "#")
}

bdsky_mcmc <- bdsky_log |>
  read_beast2_log() |>
  select(birthRate_BDSKY_Serial) |>
  mutate(
    r_naught = birthRate_BDSKY_Serial / (death_rate + sampling_rate + occurrence_rate),
    model = "bdsky"
  ) |>
  select(r_naught, model)

timtam_mcmc <- timtam_log |>
  read_beast2_log() |>
  mutate(
    r_naught = birthRate / (death_rate + sampling_rate + occurrence_rate),
    model = "timtam",
    prev_nb_mean = TimTam.prevalence.mean,
    prev_nb_size = ((TimTam.prevalence.mean ^ 2) /
                    (TimTam.prevalence.variance - TimTam.prevalence.mean))
  ) |>
  select(r_naught, model, prev_nb_mean, prev_nb_size)

## Posterior summary for the prevalence

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
  geom_histogram() +
  geom_vline(xintercept = final_prevalence)

ggsave(filename = file.path(out_dir, "prevalence-estimate.png"),
       plot = prev_gg,
       height = 10.5, width = 14.8,
       units = "cm")

## Posterior summary for R-naught

r_naught_df <- rbind(bdsky_mcmc, select(timtam_mcmc, r_naught, model))

r_naught_gg <- ggplot(
  data = r_naught_df,
  mapping = aes(x = r_naught, y = ..density..)
) +
  geom_histogram() +
  geom_vline(
    xintercept = birth_rate / (death_rate + sampling_rate + occurrence_rate)
  ) +
  facet_wrap(~model)

ggsave(filename = file.path(out_dir, "r-naught-comparison.png"),
       plot = r_naught_gg,
       height = 14.8, width = 21.0,
       units = "cm")

r_naught_summary <- r_naught_df |>
  group_by(model) |>
  summarise(
    mean = mean(r_naught),
    cred_int_lower = quantile(x = r_naught, probs = c(0.025)),
    cred_int_upper = quantile(x = r_naught, probs = c(0.975)))
write.table(x = r_naught_summary,
            file = file.path(out_dir, "r-naught-posterior-summary.csv"),
            sep = ",",
            row.names = FALSE)
