#'
#' Timtam Slam
#' ===========
#'
#' This script provides some functions that assist with setting up a
#' timtam analysis.
#'
#' Author: Alexander E. Zarebski
#'
#' Date Created: 2023-05-31
#'
#' Copyright (c) Alexander E. Zarebski, 2023
#' Email: aezarebski@gmail.com
#'
#' Provides
#' --------
#'
#' - extract_dates
#' - day_times_a
#' - fwd_day_times_a
#'
#' Dependencies
#' ------------
#'
#' - `ape` for read.FASTA

#' Extract the dates from the names of sequences in a DNAbin object
#'
#' @param seqs A DNAbin object
#'
extract_dates <- function(seqs) {
  stopifnot(class(seqs) == "DNAbin")
  strs <- names(seqs)
  r_obj <- regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", strs)
  return(regmatches(strs, r_obj))
}

#' Compute the forward times within a day to assign with sequences
#' under method A.
#'
#' @param num The number of sequences to assign times to
#'
day_times_a <- function(num) {
  method_a <- function(n, h, t) {
    tail(head(seq(from = 0, to = 1, length = n), h), t)
  }
  if (num %% 2 == 0) {
    return(method_a(num + 2, -1, -1))
  } else {
    return(method_a(num + 3, -2, -1))
  }
}

#' Compute the forward times within a day to be assigned to the
#' sequences.
#'
#' @param all_date_strs A character of sampling dates in the format
#'   "YYYY-MM-DD"
#'
fwd_day_times_a <- function(all_date_strs) {
  all_date_Date <- as.Date(all_date_strs)
  date_counts <- table(all_date_strs)
  ds_Date <- as.Date(names(date_counts))
  max_Date <- max(ds_Date)

  tmp <- round(all_date_Date - max_Date)
  tmp <- as.numeric(tmp - min(tmp))
  for (ix in seq_along(ds_Date)) {
    mask <- all_date_strs == ds_Date[ix]
    tmp[mask] <- tmp[mask] + rev(day_times_a(date_counts[ix]))
  }
  return(tmp)
}

x <- read.FASTA("sample-sequences.fasta")
all_date_strs <- extract_dates(x)
new_times <- fwd_day_times_a(all_date_strs)

## ===================================================================
## TESTS
## ===================================================================

RUN_TESTS <- TRUE

## Define a little tester function that evaluates an expression and
## stops with a suitable error message if it evaluates to false
tts_test <- function(expr) {
  if (RUN_TESTS) {
    if (!expr) {
      cat("Test expression: ", deparse(substitute(expr)), "\n")
      stop("Test failed")
    } else {
      message("Test passed")
    }
  } else {
    return(invisible())
  }
}

tts_test(day_times_a(1) == 1/3)
tts_test(all(day_times_a(2) == c(1/3, 2/3)))
tts_test(sum(abs(day_times_a(3) - c(1/5, 2/5, 3/5))) < 1e-15)
tts_test(sum(abs(day_times_a(4) - c(1/5, 2/5, 3/5, 4/5))) < 1e-15)

tts_test({
  tmp1 <- fwd_day_times_a(as.Date("2023-04-02") + c(4,3,3,0,0));
  tmp2 <- c(4 + 1/3, 3 + 2/3, 3 + 1/3, 2/3, 1/3);
  sum(abs(tmp1 - tmp2)) < 1e-15
})
