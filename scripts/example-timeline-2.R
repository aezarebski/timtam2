source("./scripts/timeline.R")


## Read out the data

timtam_xml_node <- tt_read_xml("example-2.xml")

disaster_times_bwd_rel <- tt_get_disaster_times(timtam_xml_node)
seq_times_fwd_abs <- tt_get_sequence_times(timtam_xml_node)

r0_change_times_bwd_rel <-
  tt_get_r0_change_times(timtam_xml_node)
sigma_change_times_bwd_rel <-
  tt_get_sigma_change_times(timtam_xml_node)
prop_ts_change_times_bwd_rel <-
  tt_get_prop_time_series_change_times(timtam_xml_node)
prop_seq_change_times_bwd_rel <-
  tt_get_prop_sequenced_change_times(timtam_xml_node)

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

fwd_abs_plot_xlim <- c(25.5, 31)
## fwd_abs_plot_xlim <- c(-0.1 * data_duration, 1.4 * data_duration)

ylims <- c(-1.5, 7.0)


png("example-2a.png", height = 21.0, width = 29.6, units = "cm", res = 300)

plot_background_bars(fwd_abs_plot_xlim, ylims, offset)
plot_sequence_times(seq_times_fwd_abs, num_time_labels = 2)
plot_last_sequence_indicator(seq_times_fwd_abs, vert_val = 6)
plot_time_series_times(seq_times_fwd_abs, disaster_times_bwd_rel)
plot_time_blocks(
  seq_times_fwd_abs, r0_change_times_bwd_rel,
  label_str = "R0 change times", col = "purple",
  vert_vals = c(1.0, 1.2, 1.4, 1.6)
)
plot_time_blocks(
  seq_times_fwd_abs, sigma_change_times_bwd_rel,
  label_str = "Net removal rate change times", col = "green",
  vert_vals = c(1.0, 1.2, 1.4, 1.6) + 1.0
)
plot_time_blocks(
  seq_times_fwd_abs, prop_ts_change_times_bwd_rel,
  label_str = "Time series proportion change times", col = "yellow",
  vert_vals = c(1.0, 1.2, 1.4, 1.6) + 2.0
)
plot_time_blocks(
  seq_times_fwd_abs, prop_seq_change_times_bwd_rel,
  label_str = "Time series proportion change times", col = "pink",
  vert_vals = c(1.0, 1.2, 1.4, 1.6) + 3.0
)
dev.off()
