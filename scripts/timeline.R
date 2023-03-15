library(xml2)

#' All times should be backwards and relative to the final sequence time!
#'
#' - tt_read_xml
#' - tt_get_disaster_times
#' - tt_get_sequence_times
#' - tt_get_r0_change_times
#' - plot_time_blocks
#'

plot_time_blocks <- function(fwd_abs_times, ys, label_str, ...) {
  ts <- fwd_abs_times
  if (length(ts) == 0) {
    stop("bad times given to time_blocks")
  } else if (length(ts) == 1) {
    rect(-10^10, ys[1], ts[1], ys[2], ...)
    rect(ts[1], ys[2], 10^10, ys[3], ...)
  } else {
    rect(-10^10, ys[1], ts[1], ys[2], ...)
    tmp <- TRUE
    for (ix in 1:(length(ts) - 1)) {
      if (tmp) {
        rect(ts[ix], ys[2], ts[ix + 1], ys[3], ...)
      } else {
        rect(ts[ix], ys[1], ts[ix + 1], ys[2], ...)
      }
      tmp <- !tmp
    }
    if (tmp) {
      rect(ts[ix + 1], ys[2], 10^10, ys[3], ...)
    } else {
      rect(ts[ix + 1], ys[1], 10^10, ys[2], ...)
    }
  }
  text(x = median(ts), y = ys[4], label = label_str, ...)
}

tt_read_xml <- function(file) {
  read_xml(file)
}

tt_get_disaster_times <- function(tt) {
  dts_node <-
    xml_find_first(
      x = tt,
      xpath = "//distribution[@spec=\'timtam.TimTam\']/parameter[@name='disasterTimes\']" # nolint
    )
  as.numeric(unlist(strsplit(x = xml_text(dts_node), split = " ")))
}

tt_get_origin_time <- function(tt) {
  dts_node <-
    xml_find_first(
      x = tt,
      xpath = "//distribution[@spec=\'timtam.TimTam\']/parameter[@name='originTime\']" # nolint
    )
  as.numeric(xml_text(dts_node))
}

tt_get_r0_change_times <- function(tt) {
  r0_ch_node <-
    xml_find_first(
      x = tt,
      xpath = "//distribution[@spec=\'timtam.TimTam\']/parameter[@name='r0ChangeTimes\']" # nolint
    )
  ## These times are backwards relative.
  as.numeric(unlist(strsplit(x = xml_text(r0_ch_node), split = " ")))
}

tt_get_sequence_times <- function(tt, make_bwd_rel=FALSE) {
  trt_node <-
    xml_find_first(
      x = tt,
      xpath = "//tree/trait"
    )
  ## These times are forwards and absolute.
  fwd_abs_seq_times <- as.numeric(
    sapply(
      strsplit(
        x = unlist(strsplit(x = xml_attr(x = trt_node, attr = "value"), split = ",")),
        split = "="
      ),
      function(s) {
        s[2]
      }
    )
  )
  if (!make_bwd_rel) {
    return(fwd_abs_seq_times)
  } else {
    return(max(fwd_abs_seq_times) - fwd_abs_seq_times)
  }
}


## Read out the data

disaster_times_bwd_rel <- tt_read_xml("timtam.xml") |> tt_get_disaster_times()
seq_times_fwd_abs <- tt_read_xml("timtam.xml") |> tt_get_sequence_times()
r0_change_times_bwd_rel <- tt_read_xml("timtam.xml") |> tt_get_r0_change_times()

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

## fwd_abs_plot_xlim <- c(37, 40)
fwd_abs_plot_xlim <- c(-0.1 * data_duration, 1.4 * data_duration)

ylims <- c(-1.5, 2.5)


## Plot the sequence times

data_duration <- max(seq_times_fwd_abs) - min(disaster_times_bwd_rel)

plot.new()
plot.window(xlim = fwd_abs_plot_xlim, ylim = ylims)




for (ix in seq(from = 0 - offset, to = data_duration, by = 2)) {
  rect(ix, ylims[1], ix + 1, ylims[2], col = "lightgrey", lty = 0)
}

staggered_ys <-
  head(
    rep(
      seq(from = 0, to = 0.2, length = 5),
      ceiling(length(seq_times_fwd_abs) / 5)
    ),
    length(seq_times_fwd_abs)
  )

points(seq_times_fwd_abs, staggered_ys,
  col = "darkgreen"
)

fwd_abs_unq <- unique(round(seq_times_fwd_abs, digits = 3))
tmp_mask <- round(seq(from = 1, to = length(fwd_abs_unq), length = 5))
fwd_abs_unq <- fwd_abs_unq[tmp_mask]

text(x = fwd_abs_unq, y = rep(0.4, length(fwd_abs_unq)),
     labels = as.character(fwd_abs_unq),
     col = "darkgreen")
text(x = median(fwd_abs_unq), y = 0.6,
     label = "Sequence times (forward absolute)",
     col = "darkgreen")

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
