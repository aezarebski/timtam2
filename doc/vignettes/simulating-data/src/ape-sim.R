#!/usr/bin/env Rscript
#'
VERSION <- c(0,2,1)
#' =======
#' ape-sim
#' =======
#'
#' Use the ape package to simulate the BDSCOD process from the command line.
#'
#' Simulations are conditioned on there being more than one sequenced sample. If
#' no satisfactory simulation is generated in 100 replicates then the program
#' will crash. The parameters for the simulation are read in from an XML file as
#' demonstrated below.
#'
#' Usage
#' =====
#'
#' The following will simulate without either a catastrophe or a disaster:
#'
#' $ ./ape-sim.R --verbose demo.xml
#'
#' where the \code{demo.xml} file should like the following
#'
#' <ape version="0.1.2">
#'   <configuration>
#'     <parameters birthRate="3.0"
#'                 deathRate="1.0"
#'                 samplingRate="0.5"
#'                 occurrenceRate="0.5"
#'                 duration="4.0"/>
#'     <options seed="2"
#'              writeNewick="true"
#'              makePlots="true"
#'              outputDirectory="out"
#'              simulateSequences="false"
#'              seq_agg_times=""
#'              occ_agg_times="1.0 4.0 1.0"/>
#'   </configuration>
#' </ape>
#'
#' The parameters tag can also take "nu" and "rho" as parameters if you want
#' there to be a scheduled sample at the end of the simulation. For the time
#' being the only real schema for the XML is the \code{parse_xml_configuration}
#' function.
#'
#' Help
#' ====
#'
#' $ ./ape-sim.R --help
#'
#' ChangeLog
#' =========
#'
#' - 0.2.1
#'   + Fix bug where version was not printing properly.
#'
#' - 0.2.0
#'   + Include support for time varying parameters.
#'
#' - 0.1.7
#'   + Switch to using \code{xml2} instead of \code{XML}.
#'
#' - 0.1.6
#'   + Include additional checks that configuration is valid.
#'   + Fix typo which was breaking sequence simulation.
#'
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(tidytree))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(xml2))

green_hex_colour <- "#1b9e77"
purple_hex_colour <- "#7570b3"

#' Parse a string \"FROM TO BY\" into a linear vector of values.
parse_from_to_by <- function(from_to_by_string) {
  if (from_to_by_string == "") {
    NULL
  } else {
    tmp <- as.numeric(unlist(strsplit(x = from_to_by_string, split = " ")))
    stopifnot(length(tmp) == 3)
    seq(from = tmp[1], to = tmp[2], by = tmp[3])
  }
}

#' Return a vectorised step function.
#'
#' @param ts the times at which there is a step.
#' @param vs the value of the function at the step.
#'
step_function_factory <- function(ts, vs) {
  n <- length(ts)
  if (n < 1) {
    stop("need at least one step time.")
  }
  if (length(vs) != (1 + n)) {
    stop("if there are n step times there must be n+1 values.")
  }

  step_func <- function(x) {
    m <- length(x)
    y <- numeric(m)
    mask <- logical(m)
    mask <- x < ts[1]
    y[mask] <- vs[1]
    if (n > 1) {
      for (ix in seq.int(2, n)) {
        mask <- (ts[ix-1] <= x) & (x < ts[ix])
        y[mask] <- vs[ix]
      }
    }
    mask <- x >= ts[n]
    y[mask] <- vs[n+1]
    return(y)
  }

  return(
    list(
      func = step_func,
      times = ts,
      values = vs
    )
  )
}

#' Return a vectorised step function resulting from the point-wise application
#' of a binary function to two step functions.
#'
#' @param step_func_a a step function (see \code{step_function_factory}.)
#' @param step_func_b a step function (see \code{step_function_factory}.)
#' @param binary_func a binary function.
#'
map_step_function <- function(step_func_a, step_func_b, binary_func) {
  times_c <- sort(union(step_func_a$times, step_func_b$times))
  n <- length(times_c)
  values_c <- numeric(n + 1)
  a_t_0 <- ifelse(is.finite(step_func_a$times[1]),
                  step_func_a$times[1] - 1,
                  ifelse(step_func_a$times[1] > 0,
                         0,
                         stop("invalid first time in step function a")))
  b_t_0 <- ifelse(is.finite(step_func_b$times[1]),
                  step_func_b$times[1] - 1,
                  ifelse(step_func_b$times[1] > 0,
                         0,
                         stop("invalid first time in step function b")))
  values_c[1] <- binary_func(
    step_func_a$func(a_t_0),
    step_func_b$func(b_t_0)
  )
  for (ix in seq.int(n)) {
    values_c[ix+1] <- binary_func(
      step_func_a$func(times_c[ix]),
      step_func_b$func(times_c[ix])
    )
  }
  return(step_function_factory(times_c, values_c))
}

step_func_plus <- function(sf_a, sf_b) {
  return(map_step_function(sf_a, sf_b, function(a, b) {a + b}))
}

step_func_divide <- function(sf_a, sf_b) {
  return(map_step_function(sf_a, sf_b, function(a, b) {a / b}))
}

#' Parse the rate out of the XML node
parse_rate <- function(xml_input, rate_name) {
  xml_params <- xml_find_first(xml_input, "//configuration//parameters")
  rate_attr <- xml_attr(xml_params, rate_name)
  if (grepl(pattern = "@[a-zA-Z]+", x = rate_attr)) {
    ## The rate refers to a stepFunction node so we need to extract this data
    ## and construct an appropriate step function.
    sf_node <- xml_find_first(
      xml_input,
      sprintf(
        "//stepFunction[@name='%s']",
        gsub(pattern = "@", replacement = "", x = rate_attr)
      )
    )
    ts <- as.numeric(unlist(strsplit(xml_attr(sf_node, "times"), " ")))
    vs <- as.numeric(unlist(strsplit(xml_attr(sf_node, "values"), " ")))
    return(step_function_factory(ts, vs))
  } else {
    ## To avoid having multiple representations of the rate, we can return
    ## another step function which is constant.
    return(step_function_factory(c(Inf), c(as.numeric(rate_attr), NA)))
  }
}

#' Parse the XML configuration of the simulation. See the documentation above
#' for an example.
parse_xml_configuration <- function(filepath) {
  xml_input <- read_xml(filepath)
  xml_opts <- xml_find_first(xml_input, "//configuration//options")
  xml_params <- xml_find_first(xml_input, "//configuration//parameters")

  ## We need to extract a few individual parameters so it is helpful to have a
  ## function that makes this a bit easier.
  numeric_param <- function(n) {
    tmp <- as.numeric(xml_attr(xml_params, n))
    if (is.na(tmp)) {
      stop("could not parse single parameter for ", n)
    } else {
      return(tmp)
    }
  }

  params <- list(
    birth_rate = parse_rate(xml_input, "birthRate"),           # lambda
    death_rate = parse_rate(xml_input, "deathRate"),           # mu
    sampling_rate = parse_rate(xml_input, "samplingRate"),     # psi
    occurrence_rate = parse_rate(xml_input, "occurrenceRate"), # omega
    duration = numeric_param("duration")
  )

  if (xml_has_attr(xml_params, "rho")) {
    params$rho <- numeric_param("rho")
  } else {
    params$rho <- NULL
  }

  if (xml_has_attr(xml_params, "nu")) {
    params$nu <- numeric_param("nu")
  } else {
    params$nu <- NULL
  }

  if (xml_has_attr(xml_params, "substitutionRate")) {
    params$substitution_rate <- numeric_param("substitutionRate")
  }
  if (xml_has_attr(xml_params, "seqLength")) {
    params$seq_length <- numeric_param("seqLength")
  }

  params$net_rem_rate <- step_func_plus(
    step_func_plus(
      params$death_rate,
      params$sampling_rate
    ),
    params$occurrence_rate
  )
  params$net_per_capita_event_rate <- step_func_plus(
    params$birth_rate,
    params$net_rem_rate
  )
  params$sampling_prob <- step_func_divide(
    params$sampling_rate,
    params$net_rem_rate
  )
  params$occurrence_prob <- step_func_divide(
    params$occurrence_rate,
    params$net_rem_rate
  )
  ## this is the probability that a tip that is not due to a scheduled sample
  ## will be observed (either sequenced or occurrence) rather than being a death
  ## tip.
  params$prob_observed <- step_func_plus(
    params$sampling_prob,
    params$occurrence_prob
  )
  ## this is the probability that an unscheduled sample will be sequenced.
  params$prob_sampled_given_observed <- step_func_divide(
    params$sampling_prob, params$prob_observed
  )

  opts <- list(
    seed = as.integer(xml_attr(xml_opts, "seed")),
    output_directory = xml_attr(xml_opts, "outputDirectory"),
    make_plots = as.logical(xml_attr(xml_opts, "makePlots")),
    write_newick = as.logical(xml_attr(xml_opts, "writeNewick")),
    simulate_sequences = as.logical(xml_attr(xml_opts, "simulateSequences")))

  if (xml_has_attr(xml_opts, "seq_agg_times")) {
    opts$seq_agg_times <- parse_from_to_by(xml_attr(xml_opts, "seq_agg_times"))
  } else {
    opts$seq_agg_times <- NULL
  }

  if (xml_has_attr(xml_opts, "occ_agg_times")) {
    opts$occ_agg_times <- parse_from_to_by(xml_attr(xml_opts, "occ_agg_times"))
  } else {
    opts$occ_agg_times <- NULL
  }

  list(
    params = params,
    options = opts)
}

# create parser object
parser <- ArgumentParser()

parser$add_argument(
         "--version",
         action = "store_true",
         default = FALSE,
         help = "Print version"
)
parser$add_argument(
         "-v",
         "--verbose",
         action = "store_true",
         default = FALSE,
         help = "Verbose output"
       )
parser$add_argument(
         "xml",
         type = "character",
         help = "Filepath to XML configuration"
)

#' Repeatedly run the simulation until there is a satisfactory realisation
#' obtained.
#'
#' @param params is a list of the model parameters such as birth rate.
#' @param options is a list of configuration parameters such as output
#'   directory.
#' @param is_verbose is a logical for whether to print messages.
#'
run_conditioned_simulation <- function(params, options, is_verbose) {
  max_iterations <- 100
  has_solution <- FALSE
  curr_iter <- 0
  while ((!has_solution) & (curr_iter < max_iterations)) {
    result <- tryCatch(
      run_simulation(params, options, is_verbose),
      error = function(c) "run_simulation returned a bad simulation..."
    )
    if (class(result) != "character") {
      has_solution <- TRUE
    } else {
      cat("\trepeating simulation...\n")
    }
  }
  if (options$simulate_sequences) {
    if (is_verbose) {
      cat("simulating the sequences on the reconstructed tree...\n")
    }
    result$seq_sim <- sequence_simulation(result$reconstructed_tree,
                                          params$seq_length,
                                          params$substitution_rate)
  }
  return(result)
}

#' Simulate sequences on the given tree.
#'
#' @param tree is a phylo object
#' @param len is the length of the sequence
#' @param sub_rate is the substitution rate
#'
#' @details Simulate sequences using the Jukes-Cantor model which is the default
#'   behaviour of \code{simSeq} when given a \code{phylo} object.
#'
sequence_simulation <- function(tree, len, sub_rate) {
  if (class(tree) != "phylo") {
    stop("tree must be a phylo object.")
  }
  if (is.null(len)) {
    stop("sequence length must be positive integer.")
  }
  return(simSeq(tree, l = len, rate = sub_rate))
}

#' Simulate the length of the root of the tree.
#'
#' @param params is a list of the model parameters such as birth rate.
#'
root_simulation <- function(params) {
  rate <- params$net_per_capita_event_rate$func
  times <- params$net_per_capita_event_rate$times
  if (times[1] <= 0.0) {
    stop("will not simulate root length if rate changes before or at origin.")
  }
  rl <- rexp(n = 1, rate = rate(0.0))
  ix <- 1
  while (ix <= length(times)) {
    if (rl < times[ix]) {
      return(rl)
    } else {
      rl <- times[ix] + rexp(n = 1, rate = rate(times[ix]))
      ix <- ix + 1
    }
  }
  stop("root simulation failed for an unexpected reason...")
}

#' Run a simulation and stop if something goes wrong.
#'
#' @param params is a list of the model parameters such as birth rate.
#' @param options is a list of configuration parameters such as output
#'   directory.
#' @param is_verbose is a logical for whether to print messages.
#'
run_simulation <- function(params, options, is_verbose) {
  ## Define a small amount of time that can be used to distinguish whether two
  ## numbers are equal (up to a small tolerance).
  time_eps <- 1e-10
  if (is_verbose) {
    cat("simulating tmrca and phylogeny...\n")
  }
  ## because the trees generated by \code{rlineage} start from the TMRCA rather
  ## than the origin we need to simulate the length of the root first and
  ## subtract this from the duration.
  tmrca <- root_simulation(params)
  stopifnot(params$duration > tmrca)
  ## since the birth and death rates are potentially time dependent we need to
  ## shift the time given to them to account for the delay due to the root
  ## length.
  phy <- rlineage(
    birth = function(t) {params$birth_rate$func(t + tmrca)},
    death = function(t) {params$net_rem_rate$func(t + tmrca)},
    Tmax = params$duration - tmrca,
    eps = time_eps
  )
  if (is_verbose) {
    cat("extracting tip times and labels...\n")
  }
  ## because we are going to be subsampling the tips it is useful to have some
  ## easier ways to refer to them. To simplify post-processing of the simulation
  ## data we update the tip labels to include the absolute time (forward from
  ## the origin) of the tip.
  num_tips <- length(phy$tip.label)
  tip_ix <- seq.int(num_tips)
  tip_times <- head(ape::node.depth.edgelength(phy), num_tips)
  abs_tip_times <- tip_times + tmrca
  phy$tip.label <- sprintf("%s_%f",
                           phy$tip.label,
                           abs_tip_times)
  tip_labels <- phy$tip.label
  ## TODO we can find the tips that are still extant in the simulation but it is
  ## unclear if we can use strict equality here of if we need to account for
  ## potential error in the branch lengths. It seems like this is safe...
  extant_mask <- abs_tip_times == params$duration
  extant_labels <- tip_labels[extant_mask]
  num_extant <- length(extant_labels)
  if (is_verbose) {
    cat(
      "\tthere are ",
      num_extant,
      " lineages at the end of the simulation\n"
    )
  }
  extinct_labels <- tip_labels[!extant_mask]
  num_extinct <- length(extinct_labels)
  ## We need to perform the appropriate scheduled sampling at the end of the
  ## simulation if either a rho- or nu-sample has been requested, otherwise we
  ## need to propagate the null values for the relevant variables.
  if (is.null(c(params$rho, params$nu))) {
    if (is_verbose) {
      cat("skipping rho-sampling\n")
      cat("skipping nu-sampling\n")
    }
    num_rho_sampled <- NULL
    rho_sampled_labels <- NULL
    num_nu_sampled <- NULL
    nu_sampled_labels <- NULL
  } else if (is.null(params$nu)) {
    if (is_verbose) {
      cat("performing rho sampling...\n")
      cat("skipping nu-sampling\n")
    }
    num_rho_sampled <- rbinom(
      n = 1,
      size = num_extant,
      prob = params$rho
    )
    rho_sampled_labels <- sample(
      x = extant_labels,
      size = num_rho_sampled,
      replace = FALSE
    )
    num_nu_sampled <- NULL
    nu_sampled_labels <- NULL
  } else {
    if (is_verbose) {
      cat("skipping rho-sampling\n")
      cat("performing nu sampling...\n")
    }
    num_rho_sampled <- NULL
    rho_sampled_labels <- NULL
    num_nu_sampled <- rbinom(
      n = 1,
      size = num_extant,
      prob = params$nu
    )
    nu_sampled_labels <- sample(
      x = extant_labels,
      size = num_nu_sampled,
      replace = FALSE
    )
  }

  ## Pre-allocate the logical masks (all initially \code{FALSE}) to hold the
  ## outcome of each extinct lineage. The \code{outcome} vector is used as a
  ## reference for the fate of each tip in the tree.
  observed_mask <- rep(FALSE, num_tips)
  unsched_seq_mask <- rep(FALSE, num_tips)
  unsched_unseq_mask <- rep(FALSE, num_tips)
  outcome <- rep("", num_tips)

  for (ix in seq.int(num_tips)) {
    if (!extant_mask[ix]) {
      ## lineage ended before the present so must be "sampling", "occurrence" or
      ## "death".
      is_observed <- rbinom(
        n = 1,
        size = 1,
        prob = params$prob_observed$func(abs_tip_times[ix])
      )
      if (is_observed == 1) {
        observed_mask[ix] <- TRUE
        is_seq <- rbinom(
          n = 1,
          size = 1,
          prob = params$prob_sampled_given_observed$func(abs_tip_times[ix])
        )
        if (is_seq == 1) {
          outcome[ix] <- "sampling"
          unsched_seq_mask[ix] <- TRUE
        } else {
          outcome[ix] <- "occurrence"
          unsched_unseq_mask[ix] <- TRUE
        }
      } else {
        outcome[ix] <- "death"
      }
    } else {
      ## lineage persisted to the present so must be "rho", "nu", or "extant".
      if (is.null(c(rho_sampled_labels, nu_sampled_labels))) {
        outcome[ix] <- "extant"
      } else {
        tl <- tip_labels[ix]
        if (is.element(tl, rho_sampled_labels)) {
          outcome[ix] <- "rho"
        } else if (is.element(tl, nu_sampled_labels)) {
          outcome[ix] <- "nu"
        } else {
          outcome[ix] <- "extant"
        }
        rm(tl)
      }
    }
  }

  num_observed <- sum(observed_mask)
  num_sampled <- sum(unsched_seq_mask)
  observed_labels <- tip_labels[observed_mask]
  sampling_labels <- tip_labels[unsched_seq_mask]
  occurrence_labels <- tip_labels[unsched_unseq_mask]
  num_occurrences <- length(occurrence_labels)
  unobserved_labels <- setdiff(extinct_labels, observed_labels)

  ## We then drop the lineages that did not result in a sampled/rho lineage and
  ## the remaining tree is the reconstructed tree.
  rt_tip_labels <- c(sampling_labels, rho_sampled_labels)
  num_rt_tips <- length(rt_tip_labels)
  if (num_rt_tips < 2) {
    stop(paste0(c("\n\n",
                  rep("-", 60),
                  "\n\n\tSIMULATION HAS LESS THAN TWO SEQUENCED SAMPLES!\n\n",
                  rep("-", 60),
                  "\n"),
                collapse = ""))
  }
  rt <- ape::keep.tip(phy, rt_tip_labels)
  rt_tip_and_node_depths <- ape::node.depth.edgelength(rt)
  rt_tip_depths <- head(rt_tip_and_node_depths, num_rt_tips)
  rt_node_depths <- tail(rt_tip_and_node_depths, -num_rt_tips)
  ## The depths of the nodes and tips is relative to the TMRCA so we need to
  ## adjust for this when computing the exact times that they occurred at. Since
  ## the TMRCA in the reconstructed tree may be different to the one of the
  ## whole tree we cannot simply re-use the previous rootlength.
  true_first_sample_time <-
    min(tip_times[outcome == "sampling" | outcome == "rho"])
  depth_first_in_rt <- min(rt_tip_depths)
  rt_offset <- true_first_sample_time - depth_first_in_rt
  ## we can put all of this information into a single dataframe which is more
  ## convenient for export and subsequent usage. Note that the times here are
  ## relative to the TMRCA of the whole transmission tree, not the reconstructed
  ## tree or the origin. This is why we use the tip times and not the \code{rt_}
  ## times.
  event_times_df <- data.frame(
    time = c(-tmrca,
             tip_times[outcome == "occurrence"],
             tip_times[outcome == "sampling"],
             tip_times[outcome == "rho"],
             tip_times[outcome == "nu"],
             rt_node_depths + rt_offset
             ),
    event = c("origin",
              rep(
                c("occurrence", "sampling", "rho", "nu", "birth"),
                times = c(num_occurrences,
                          num_sampled,
                          ifelse(is.null(num_rho_sampled), 0, num_rho_sampled),
                          ifelse(is.null(num_nu_sampled), 0, num_nu_sampled),
                          ifelse(is.null(num_rho_sampled), num_sampled - 1, num_rt_tips - 1)))))

  ## Compute the prevalence of infection at the end of the simulation. If there
  ## is no rho-sampling and no nu-sampling then this is just the number of
  ## extant lineages. Otherwise, we need to subtract the number that were
  ## rho-sampled or nu-sampled, depending on which was carried out, and remove
  ## this from the number of extant lineages.
  if (is.null(c(params$rho, params$nu))) {
    final_prevalence <- num_extant
  } else if (is.null(params$nu)) {
    final_prevalence <- num_extant - num_rho_sampled
  } else {
    final_prevalence <- num_extant - num_nu_sampled
  }
  if (is_verbose) {
    cat("checking output from run_simulation...\n")
  }

  ## TODO The aggregation code below should probably be refactored and applied
  ## as a post-simulation step rather than being included here.

  ## If there is aggregation that needs to be carried out then this needs to be
  ## done now just before the simulation results are returned. This operation is
  ## only defined when there are no rho-samples and no nu-samples at the end of
  ## the simulation (which should be enforced by the validation of the
  ## parameters.) Here we construct the \code{agg_event_times_df} which contains
  ## the data after accounting for the aggregation.
  if (is.null(c(params$rho, params$nu)) &&
      (!is.null(options$occ_agg_times) | !is.null(options$seq_agg_times))) {
    ## There is a little extra book keeping involved because the times are
    ## relative to the TMRCA of the whole transmission tree rather than the
    ## origin. To make sense of the requested aggregation times we need to
    ## adjust them to be relative to the TMRCA of the whole transmission tree.
    ## This is why we subtract the TMRCA from the given aggregation times. We
    ## add a column for size to count the number of individuals that were
    ## removed in the mock scheduled event (resulting from the aggregation).
    origin_event_row <- filter(event_times_df, event == "origin")
    origin_event_row$size <- NA
    birth_rows <- filter(event_times_df, event == "birth")
    birth_rows$size <- NA

    if (!is.null(options$seq_agg_times)) {
      ## there are sequence aggregation times so we need to extract these rows
      ## in a suitable data frame. Recall that the times in
      ## \code{event_times_df} are relative to the TMRCA of the whole
      ## transmission tree. We construct a data frame, \code{seq_rows} to
      ## represent the aggregated counts and their times.
      sampling_times <- select(filter(event_times_df, event == "sampling"), time)
      tmrca_rel_seq_agg_times <- options$seq_agg_times - tmrca
      num_agg_obs <- length(tmrca_rel_seq_agg_times) - 1
      if (max(options$seq_agg_times) < max(sampling_times)) {
        stop(
          "There are sequenced samples after the last given aggregation time.",
          "It is unclear how to account for the births related to these",
          "sequences so you probably want to include another aggregation point",
          "to capture them."
        )
      }
      tmp_time <- numeric(num_agg_obs)
      tmp_size <- numeric(num_agg_obs)
      for (ix in seq.int(num_agg_obs)) {
        tmp_time[ix] <- tmrca_rel_seq_agg_times[ix+1]
        tmp_size[ix] <- sampling_times |>
          filter(tmrca_rel_seq_agg_times[ix] < time,
                 time <= tmrca_rel_seq_agg_times[ix+1]) |>
          nrow()
      }
      seq_rows <- data.frame(time = tmp_time,
                             event = "rho",
                             size = tmp_size)
    } else {
      ## If there are not any sequence aggregation times then the sequenced
      ## samples still need to be represented by sampling events so we pull
      ## these values out and store them in an appropriate data frame for
      ## subsequent use.
      seq_rows <- event_times_df[event_times_df$event == "sampling", ]
      seq_rows$size <- NA
    }

    if (!is.null(options$occ_agg_times)) {
      occurrence_times <- select(filter(event_times_df, event == "occurrence"), time)
      tmrca_rel_occ_agg_times <- options$occ_agg_times - tmrca
      num_agg_obs <- length(tmrca_rel_occ_agg_times) - 1
      tmp_time <- numeric(num_agg_obs)
      tmp_size <- numeric(num_agg_obs)
      for (ix in seq.int(num_agg_obs)) {
        tmp_time[ix] <- tmrca_rel_occ_agg_times[ix+1]
        tmp_size[ix] <- occurrence_times |>
          filter(tmrca_rel_occ_agg_times[ix] < time,
                 time <= tmrca_rel_occ_agg_times[ix+1]) |>
          nrow()
      }
      unseq_rows <- data.frame(time = tmp_time,
                               event = "nu",
                               size = tmp_size)
    } else {
      unseq_rows <- event_times_df[event_times_df$event == "occurrence", ]
      unseq_rows$size <- NA
    }
    agg_event_times_df <- rbind(origin_event_row,
                                unseq_rows,
                                seq_rows,
                                birth_rows)
  } else {
    ## no aggregation has been called for so this variable gets set to
    ## \code{NULL}
    agg_event_times_df <- NULL
  }

  return(list(
    event_times_df = event_times_df,
    aggregated_event_times_df = agg_event_times_df,
    final_prevalence = final_prevalence,
    phylo = phy,
    reconstructed_tree = rt,
    outcome = outcome,
    tip_ix = tip_ix,
    num_extinct = num_extinct,
    num_extant = num_extant,
    num_sampled = num_sampled,
    num_rho_sampled = num_rho_sampled,
    num_nu_sampled = num_nu_sampled,
    num_observed = num_observed,
    num_occurrences = num_occurrences))
}

#' Make some nice plots and write them to the output directory.
write_plot <- function(simulation_results,
                       parameters,
                       output_directory,
                       is_verbose) {

  is_observed <- is.element(
    simulation_results$outcome,
    c("sampling", "occurrence", "rho", "nu")
  )

  tip_annotations <- tibble(
    node = simulation_results$tip_ix,
    outcome = simulation_results$outcome,
    is_observed = is_observed
  )

  tr <- tidytree::treedata(
                    phylo = simulation_results$phylo,
                    data = tip_annotations
                  )

  colourful_tree_plot <- ggplot(tr, mapping = aes(x, y)) +
    geom_tippoint(mapping = aes(colour = outcome, shape = is_observed),
                  size = 3) +
    geom_nodepoint() +
    geom_vline(data = simulation_results$event_times_df,
               aes(xintercept = time, colour = event)) +
    geom_tree() +
    labs(colour = "Tip outcome",
         shape = "Observed") +
    theme_tree2(legend.position = "top")

  neat_tree_plot <- ggplot(tr, mapping = aes(x, y)) +
    geom_tree(alpha = 0.3) +
    geom_tippoint(mapping = aes(colour = is_observed),
                  size = 2.5) +
    scale_colour_manual(
      values = c("#bdbdbd", green_hex_colour),
      labels = c("TRUE" = "Observed", "FALSE" = "Not observed")) +
    labs(colour = NULL) +
    theme_tree(legend.position = "top")

  if (is_verbose) {
    cat("writing visualistion to file...\n")
  }
  ggsave(
    filename = paste(
      output_directory,
      "ape-simulation-figure.png",
      sep = "/"
    ),
    plot = colourful_tree_plot,
    width = 25,
    height = 15,
    units = "cm"
  )
  saveRDS(object = list(colourful_tree_plot, neat_tree_plot),
          file = paste(
            output_directory,
            "ape-sim-figures-1.rds",
            sep = "/")
          )
}

write_aggregated_plot <- function(simulation_results,
                                  parameters,
                                  output_directory,
                                  is_verbose) {
  plot_df <- filter(
    simulation_results$aggregated_event_times_df,
    is.element(event, c("rho", "nu"))
  )
  g <- ggplot(
    data = plot_df,
    mapping = aes(x = time, y = size, colour = event)
  ) +
    geom_line() +
    geom_point() +
    scale_colour_manual(
      breaks = c("nu", "rho"),
      values = c(purple_hex_colour, green_hex_colour)
    ) +
    labs(y = "Count", x = "Day", colour = "Observation type") +
    theme_classic()
  ggsave(
    filename = paste(
      output_directory,
      "ape-simulation-figure-aggregated.png",
      sep = "/"
    ),
    plot = g,
    width = 14.8,
    height = 10.5,
    units = "cm"
  )
  g2 <- ggplot(
    data = plot_df,
    mapping = aes(x = time, y = size, shape = event)
  ) +
    geom_line(
      colour = purple_hex_colour
    ) +
    geom_point(
      colour = purple_hex_colour,
      size = 3
    ) +
    scale_shape_manual(
      values = c(1, 2),
      labels = c("nu" = "Weekly unsequenced", "rho" = "Daily sequenced")
    ) +
    scale_y_sqrt() +
    labs(y = "Count (square root scale)", x = "Day", shape = "Observation type") +
    theme_classic() +
    theme(legend.position = c(0.2, 0.8))
  saveRDS(object = list(g2),
          file = paste(
            output_directory,
            "ape-sim-figures-2.rds",
            sep = "/")
          )
}

#' Predicate for the configuration read from XML being valid.
configuration_is_valid <- function(config) {
  params <- config$params
  opts <- config$options
  if ((!is.null(params$rho)) && (!is.null(params$seq_agg_times))) {
    return(FALSE)
  }
  if (is.element("rho", names(params))) {
    if (params$rho > 1.0 || 0.0 > params$rho) {
      warning("rho parameter should be a probability")
      return(FALSE)
    }
    if (!is.null(opts$seq_agg_times)) {
      warning("cannot have sequence aggregation times and rho-sampling.")
      return(FALSE)
    }
  }
  if (is.element("nu", names(params))) {
    if (params$nu > 1.0 || 0.0 > params$nu) {
      warning("nu parameter should be a probability")
      return(FALSE)
    }
    if (!is.null(opts$occ_agg_times)) {
      warning("cannot have unsequenced aggregation times and nu-sampling.")
      return(FALSE)
    }
  }
  ## If we are being asked to use a particular output directory it should exist.
  if (!dir.exists(opts$output_directory)) {
    warning("missing output directory: ", opts$output_directory)
    return(FALSE)
  }
  ## If we are simulating sequences we need the data to do this.
  if (opts$simulate_sequences) {
    if (is.null(params$substitution_rate)) {
      stop("cannot simulate sequences without a substitution rate.")
    }
    if (is.null(params$seq_length)) {
      stop("cannot simulate sequences without a sequence length.")
    }
  }
  return(TRUE)
}

main <- function(args, config) {
  if (args$version) {
    ## TODO when requesting the version number you also need to provide a valid
    ## XML which is silly. Because \code{main} takes both the args and the
    ## config this is an annoying fix so I'm leaving it for now.
    cat(paste(c("ape-sim-",
                paste(as.character(VERSION), collapse = "."),
                "\n"), collapse = ""))
    return(0)
  }

  if (!configuration_is_valid(config)) {
    print(config)
    stop("THE CONFIGURATION IS NOT VALID.")
  }

  if (args$verbose) {
    cat("reading parameters from", args$xml, "\n")
  }
  ## We set the seed of the PRNG so that each time this program is called the
  ## results can be reproduced exactly.
  set.seed(config$options$seed)

  ## Awkward because you cannot assign NULL via ifelse so you need to have a
  ## full conditional to do this.
  if (is.null(config$params$rho)) {
    maybe_rho <-  NULL
  } else {
    maybe_rho <- config$params$rho
  }
  if (is.null(config$params$nu)) {
    maybe_nu <- NULL
  } else {
    maybe_nu <- config$params$rho
  }

  if (args$verbose) {
    if (is.null(maybe_rho)) {
      cat("without rho sampling at the end of the simulation\n")
    } else {
      cat("with rho sampling at the end of the simulation\n")
    }
    if (is.null(maybe_nu)) {
      cat("without nu sampling at the end of the simulation\n")
    } else {
      cat("with nu sampling at the end of the simulation\n")
    }
  }

  ## the default simulator can produce simulations where the tree is not valid
  ## (because it only has a single tip) so we use the conditioned simulator to
  ## repeat the simulation if this happens.
  sim_result <- run_conditioned_simulation(config$params, config$options, args$verbose)
  if (file.access(config$options$output_directory, mode = 2) != 0) {
    stop("Cannot write to output directory: ", config$options$output_directory)
  } else {
    output_filepath <- function(n) {
      paste(config$options$output_directory, n, sep = "/")
    }
    ## record the events that happened in the simulation
    if (args$verbose) {
      cat("writing output to csv...\n")
    }
    write.table(x = sim_result$event_times_df,
                file = output_filepath("ape-sim-event-times.csv"),
                sep = ",",
                row.names = FALSE)

    if (config$options$write_newick) {
      if (args$verbose) {
        cat("writing the newick for the reconstructed tree...\n")
      }
      write.tree(phy = sim_result$reconstructed_tree,
                 file = output_filepath("ape-sim-reconstructed-tree.newick"))
    }
    if (config$options$simulate_sequences) {
      if (args$verbose) {
        cat("writing the sequence simulation...\n")
      }
      write.phyDat(sim_result$seq_sim,
                   file = output_filepath("ape-sim-sequences.fasta"),
                   format = "fasta")
    }
    ## if there are aggregated simulation results then these should also be
    ## recorded.
    if (!is.null(sim_result$aggregated_event_times_df)) {
      if (args$verbose) {
        cat("writing aggregated output to csv...\n")
      }
      write.table(x = sim_result$aggregated_event_times_df,
                  file = output_filepath("ape-sim-aggregated-event-times.csv"),
                  sep = ",",
                  row.names = FALSE)
    }
    ## record the final prevalence accounting for the number of lineages that
    ## may have been removed in a rho sample at the very end of the simulation.
    if (args$verbose) {
      cat("writing final prevalence to json...\n")
    }
    jsonlite::write_json(
                x = sim_result$final_prevalence,
                path = output_filepath("ape-sim-final-prevalence.json")
              )
  }

  if (config$options$make_plots) {
    write_plot(sim_result, config$params, config$options$output_directory, args$verbose)
    if (!is.null(sim_result$aggregated_event_times_df)) {
      write_aggregated_plot(sim_result, config$params, config$options$output_directory, args$verbose)
    }
  }
}

if (!interactive()) {
  args <- parser$parse_args()
  config <- parse_xml_configuration(args$xml)
  main(args, config)
} else {
  demo_xml_1 <- "<ape version=\"0.1.2\"><configuration><parameters birthRate=\"3.0\" deathRate=\"1.0\" samplingRate=\"0.5\" occurrenceRate=\"0.5\" duration=\"4.0\" /><options seed=\"2\" writeNewick=\"true\" makePlots=\"true\" outputDirectory=\"out\" simulateSequences=\"false\" seq_agg_times=\"\" occ_agg_times=\"1.0 4.0 1.0\" /></configuration></ape>"
  demo_xml_2 <- "<ape version=\"0.2.0\"><stepFunction name=\"stepBirthRate\" times=\"2.0 3.0\" values=\"3.0 2.5 2.0\" /><configuration><parameters birthRate=\"@stepBirthRate\" deathRate=\"1.0\" samplingRate=\"0.5\" occurrenceRate=\"0.5\" duration=\"4.0\" /><options seed=\"2\" writeNewick=\"true\" makePlots=\"true\" outputDirectory=\"out\" simulateSequences=\"false\" seq_agg_times=\"\" occ_agg_times=\"1.0 4.0 1.0\" /></configuration></ape>"
  demo_xml_3 <- "<ape version=\"0.1.2\"><configuration><parameters birthRate=\"3.0\" deathRate=\"1.0\" samplingRate=\"0.5\" occurrenceRate=\"0.5\" duration=\"4.0\" rho=\"0.5\" /><options seed=\"2\" writeNewick=\"true\" makePlots=\"true\" outputDirectory=\"out\" simulateSequences=\"false\" /></configuration></ape>"
  tmp_file <- tempfile()
  writeLines(
    text = demo_xml_3,
    con = tmp_file
  )
  args <- list(
    version = FALSE,
    verbose = TRUE,
    xml = tmp_file
  )
  config <- parse_xml_configuration(args$xml)
  main(args, config)

}

# LocalWords:  TMRCA
