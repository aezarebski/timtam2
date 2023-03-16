library(xml2)

#' All times should be backwards and relative to the final sequence time!
#'
#' - tt_read_xml
#' - tt_get_disaster_times
#' - tt_get_sequence_times
#' - tt_get_r0_change_times
#'
#' - plot_time_blocks
#' - plot_background_bars
#' - plot_sequence_times
#'

#' Plot the sequence times
#'
#' @param fwd_abs_seq_times the forwards absolute times of sequences.
#' @param vert_vals the vertical position of the elements of the plot.
#' @param num_time_labels the number of times to label.
#'
plot_sequence_times <- function(fwd_abs_seq_times,
                                vert_vals = c(0, 0.2, 0.4, 0.6),
                                num_time_labels = 5) {
  staggered_ys <-
    head(
      rep(
        seq(from = vert_vals[1], to = vert_vals[2], length = 5),
        ceiling(length(fwd_abs_seq_times) / 5)
      ),
      length(fwd_abs_seq_times)
    )

  points(fwd_abs_seq_times, staggered_ys,
    col = "darkgreen"
  )

  fwd_abs_unq <- unique(round(fwd_abs_seq_times, digits = 3))
  tmp_mask <- round(seq(
    from = 1, to = length(fwd_abs_unq),
    length = num_time_labels
  ))
  fwd_abs_unq <- fwd_abs_unq[tmp_mask]

  text(
    x = fwd_abs_unq, y = rep(vert_vals[3], length(fwd_abs_unq)),
    labels = as.character(fwd_abs_unq),
    col = "darkgreen"
  )
  text(
    x = median(fwd_abs_unq), y = vert_vals[4],
    label = "Sequence times (forward absolute)",
    col = "darkgreen"
  )
}

#' Create plot frame and add grey bars to indicate units of time.
#'
#' @param fwd_abs_xlims the interval of forwards absolute times to plot
#' @param ylims the vertical limits for the figure
#' @param offset the offset between the time of the last sequence and the end of
#'   that day.
#'
plot_background_bars <- function(fwd_abs_xlims, ylims, offset) {
  plot.new()
  plot.window(xlim = fwd_abs_xlims, ylim = ylims)

  for (ix in seq(from = 0 - offset, to = max(fwd_abs_xlims), by = 2)) {
    rect(ix, ylims[1], ix + 1, ylims[2], col = "lightgrey", lty = 0)
  }
}

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
