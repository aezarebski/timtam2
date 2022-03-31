library(dplyr)
library(reshape2)
library(ggplot2)
library(magrittr)
library(purrr)

out_dir <- "out"
bdsky_log <- file.path(out_dir, "ex2-bdsky-serial.log")
birth_rate <- NULL # 4.5 then 2.5 changing at time 2.
birth_rate_change_times <- NULL
death_rate <- 1.0
sampling_rate <- 0.5
occurrence_rate <- 0.5
final_prevalence <- 749


read_beast2_log <- function(filename) {
  read.csv(filename, sep = "\t", comment.char = "#")
}


bdsky_mcmc <- bdsky_log |>
  read_beast2_log() |>
  select(matches("reproductiveNumber*"))


## TODO Obviously everything below here needs to be fixed but this demonstrates
## the type of figure that should be created.


tmp_bdsky <- data.frame(
  t = c(0, 2.5-1e-6, 2.5 + 1e-6, 5.0),
  y = c(1.739, 1.739, 1.205, 1.205),
  ymin = c(1.558114, 1.558114, 1.134191, 1.134191),
  ymax = c(1.914309, 1.914309, 1.275801, 1.275801)
)

tmp_sim <- data.frame(
  t = c(0, 2.0 - 1e-6, 2.0 + 1e-6, 5.0),
  y = c(4.5 / 2.0, 4.5 / 2.0, 2.5 / 2.0, 2.5 / 2.0)
)

tmp_gg <- ggplot() +
  geom_ribbon(
    data = tmp_bdsky,
    mapping = aes(x = t, ymin = ymin, ymax = ymax),
    alpha = 0.1
  ) +
  geom_line(
    data = tmp_bdsky,
    mapping = aes(x = t, y = y),
  ) +
  geom_line(
    data = tmp_sim,
    mapping = aes(x = t, y = y),
    linetype = "dashed"
  )

ggsave(filename = file.path(out_dir, "demonstration-plot.png"),
       plot = tmp_gg,
       height = 10.5, width = 14.8,
       ## A5 height = 14.8, width = 21.0,
       ## A6 height = 10.5, width = 14.8,
       ## A7 height = 7.4, width = 10.5,
       units = "cm")
