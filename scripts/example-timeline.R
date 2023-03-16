source("./scripts/timeline.R")


## Read out the data

timtam_xml_node <- tt_read_xml("timtam.xml")

disaster_times_bwd_rel <- tt_get_disaster_times(timtam_xml_node)
seq_times_fwd_abs <- tt_get_sequence_times(timtam_xml_node)
r0_change_times_bwd_rel <- tt_get_r0_change_times(timtam_xml_node)

#'
#' this needs to be calculated by hand depending upon what time of day you are
#' assuming the last sequence was collected on!!!
#'
#' This specifies the amount of units of time (days) between the time of the
#' last sequence and the end of that day.
#'
offset <- 0.0 # calculate by hand!

#' Specify the range of times we are going to care about. This is how you
#' control zooming into small regions of time to understand exactly when events
#' are occuring.
#'

data_duration <- max(seq_times_fwd_abs) - min(disaster_times_bwd_rel)

## fwd_abs_plot_xlim <- c(37, 40)
fwd_abs_plot_xlim <- c(-0.1 * data_duration, 1.4 * data_duration)

ylims <- c(-1.5, 2.5)

## Create the canvas to draw the data on.

plot_background_bars(fwd_abs_plot_xlim, ylims, offset)

## Plot the sequence times

plot_sequence_times(seq_times_fwd_abs)

## Add an indication of the relative zero time.

fwd_abs_zero <- max(seq_times_fwd_abs)
abline(v = fwd_abs_zero, lty = "dashed")
text(
  x = fwd_abs_zero, y = 2,
  label = sprintf("Last sequence:\n\t\t\t%.3f (forward absolute)\n\t\t\t%.3f (backwards relative)", fwd_abs_zero, 0) # nolint
)

## Draw the time series of cases

points(x = fwd_abs_zero - disaster_times_bwd_rel, y = rep(-0.5, length(fwd_abs_zero - disaster_times_bwd_rel)), col = "darkred")

bwd_rel_mask <- round(seq(from = 1, to = length(disaster_times_bwd_rel), length = 5))
bwd_rel_times <- disaster_times_bwd_rel[bwd_rel_mask]
text(
  x = fwd_abs_zero - bwd_rel_times, y = -0.8,
  labels = as.character(round(bwd_rel_times, digits = 3)),
  col = "darkred"
)
text(x = median(fwd_abs_zero - bwd_rel_times), y = -1.0,
     labels = "Time series times (backwards relative)",
     col = "darkred")

## Draw the times of any changes to R0

plot_time_blocks(fwd_abs_zero - r0_change_times_bwd_rel, c(1.0, 1.2, 1.4, 1.6), "R0 change times", col = "purple")
