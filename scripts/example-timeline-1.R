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

plot_last_sequence_indicator(seq_times_fwd_abs)

## Draw the time series of cases

plot_time_series_times(seq_times_fwd_abs, disaster_times_bwd_rel)

## Draw the times of any changes to R0

plot_time_blocks(seq_times_fwd_abs, r0_change_times_bwd_rel,
                 label_str = "R0 change times", col = "purple")
