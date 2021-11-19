#!/usr/bin/env Rscript
#'
VERSION <- c(0,1,2)
#' ape-sim-0.1.2
#' =============
#'
#' Use the ape package to simulate the BDSCOD process from the command line.
#'
#' Simulations are conditioned on there being more than one sequenced sample. If
#' no satisfactory simulation is generated in 100 replicates then the program
#' will crash. By default this avoids a population sample at the end of the
#' simulation but there is a command line argument if you want to include this.
#'
#' Usage
#' -----
#'
#' The following will simulate without either a catastrophe or a disaster:
#'
#' $ ./ape-sim.R --seed 1 -p my-params.json -o out --duration 3.0
#'
#' where the \code{my-params.json} file should like the following
#'
#' {
#'   "birthRate": 3.0,
#'   "deathRate": 1.0,
#'   "samplingRate": 0.5,
#'   "occurrenceRate": 0.5
#' }
#'
#' If you want to include a rho-sample at the end there is a --rho flag for
#' that. Or if you want to generate a figure showing the simulated data append
#' the --make-plots flag.
#'
#' Help
#' ----
#'
#' $ ./ape-sim.R --help
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

green_hex_colour <- "#1b9e77"
purple_hex_colour <- "#7570b3"

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
         "-s",
         "--seed",
         type = "integer",
         default = 1,
         help = "PRNG seed"
       )
parser$add_argument(
         "-p",
         "--parameters",
         type = "character",
         help = "Filepath to parameters JSON"
       )
parser$add_argument(
         "-d",
         "--duration",
         type = "double",
         help = "Simulation duration"
       )
parser$add_argument(
         "-o",
         "--output-directory",
         type = "character",
         help = "Path to write output to"
       )
parser$add_argument(
         "--make-plots",
         action = "store_true",
         default = FALSE,
         help = "Generate plots"
       )
parser$add_argument(
         "--write-newick",
         action = "store_true",
         default = FALSE,
         help = "Write a Newick string representation of the reconstructed tree to file."
       )
parser$add_argument(
         "-r",
         "--rho",
         type = "double",
         default = -1.0,
         help = "Scheduled sample probability at the end of the simulation.")
parser$add_argument(
         "--seq-agg-times",
         type = "character",
         default = "",
         help = "Specification of aggregation times for sequenced samples: \"FROM TO BY\". These values get read as three numbers then form the arguments for the seq function.")
parser$add_argument(
         "--occ-agg-times",
         type = "character",
         default = "",
         help = "Specification of aggregation times for unsequenced samples (occurrence data). See the details of --seq-agg-times.")
parser$add_argument(
         "--simulate-sequences",
         action = "store_true",
         default = FALSE,
         help = "Simulate sequences for each of the leaves of the reconstructed tree. This makes use of the simSeq function from phangorn. If this parameter is given, then it is assumed that there will be a substitutionRate given in the parameters JSON.")

read_parameters <- function(parameter_filepath,
                            sim_duration,
                            maybe_rho,
                            maybe_seq_from_to_by,
                            maybe_occ_from_to_by,
                            simulate_sequences,
                            is_verbose) {
  if (!file.exists(parameter_filepath)) {
    stop("Cannot find parameter file: ", parameter_filepath)
  } else if (sim_duration <= 0.0) {
    stop("Need a positive duration")
  } else if (!(is.null(maybe_rho) | is.numeric(maybe_rho))) {
    stop("Need a maybe numeric value for maybe_rho")
  } else {
    params <- jsonlite::read_json(parameter_filepath)
    ## it will be useful to have some other ways to talk about the parameters so
    ## we compute a couple of other views.
    params$net_rem_rate <-
      params$deathRate + params$samplingRate + params$occurrenceRate
    params$net_per_capita_event_rate <- params$birthRate + params$net_rem_rate
    params$sampling_prob <- params$samplingRate / params$net_rem_rate
    params$occurrence_prob <- params$occurrenceRate / params$net_rem_rate
    params$prob_observed <- params$sampling_prob + params$occurrence_prob
    params$prob_sampled_given_observed <-
      params$sampling_prob / params$prob_observed
    params$duration <- sim_duration
    ## we need to handle the possibility that there is a valid rho.
    if (is.null(maybe_rho)) {
      params$rho <- maybe_rho
    } else if (0 < maybe_rho && maybe_rho < 1) {
      params$rho <- maybe_rho
    } else {
      stop("invalid rho argument: ", maybe_rho)
    }
    ## we need to parse the aggregation times if they have been given.
    if (is.null(maybe_seq_from_to_by)) {
      params$seq_agg_times <- NULL
    } else {
      params$seq_agg_times <- parse_from_to_by(maybe_seq_from_to_by)
    }
    if (is.null(maybe_occ_from_to_by)) {
      params$occ_agg_times <- NULL
    } else {
      params$occ_agg_times <- parse_from_to_by(maybe_occ_from_to_by)
    }
    params$simulate_sequences <- simulate_sequences
    return(params)
  }
}

run_conditioned_simulation <- function(params, is_verbose) {
  max_iterations <- 100
  has_solution <- FALSE
  curr_iter <- 0
  while ((!has_solution) & (curr_iter < max_iterations)) {
    result <- tryCatch(
      run_simulation(params, is_verbose),
      error = function(c) "run_simulation returned a bad simulation..."
    )
    if (class(result) != "character") {
      has_solution <- TRUE
    } else {
      cat("\trepeating simulation...\n")
    }
  }
  if (params$simulate_sequences) {
    if (is_verbose) {
      cat("simulating the sequences on the reconstructed tree...\n")
    }
    result$seq_sim <- sequence_simulation(result$reconstructed_tree,
                                          params$substitutionRate)
  }
  return(result)
}

sequence_simulation <- function(tr, sub_rate) {
  return(simSeq(tr, rate = sub_rate))
}

run_simulation <- function(params, is_verbose) {
  time_eps <- 1e-6
  if (is_verbose) {
    cat("simulating tmrca and phylogeny...\n")
  }
  ## because the trees generated by rlineage start from the TMRCA rather than
  ## the origin we need to simulate the length of the root first and subtract
  ## this from the duration.
  tmrca <- rexp(n = 1, rate = params$net_per_capita_event_rate)
  stopifnot(params$duration > tmrca)
  phy <- rlineage(
    params$birthRate,
    params$net_rem_rate,
    Tmax = params$duration - tmrca,
    eps = time_eps
  )
  if (is_verbose) {
    cat("extracting tip times and labels...\n")
  }
  ## because we are going to be subsampling the tips it is useful to have some
  ## easier ways to refer to them. To make simplify post-processing of the
  ## simulation data we update the tip labels to include the absolute time
  ## (forward from the origin) of the tip.
  num_tips <- length(phy$tip.label)
  tip_ix <- seq.int(num_tips)
  tip_times <- head(ape::node.depth.edgelength(phy), num_tips)
  phy$tip.label <- sprintf("%s_%f",
                           phy$tip.label,
                           tip_times + tmrca)
  tip_labels <- phy$tip.label
  ## TODO we can find the tips that are still extant in the simulation but it is
  ## unclear if we can use strict equality here of if we need to account for
  ## potential error in the branch lengths. It seems like this is safe...
  extant_mask <- tip_times + tmrca == params$duration
  extant_labels <- tip_labels[extant_mask]
  num_extant <- length(extant_labels)
  if (is_verbose) {
    cat("\tthere appear to be ", num_extant, " lineages at the end of the simulation\n")
  }
  ## We need to do a rho sample at the end of the duration if one has been
  ## requested, otherwise we need to propagate the null values.
  if (is.null(params$rho)) {
    if (is_verbose) {
      cat("skipping rho sampling\n")
    }
    num_rho_sampled <- NULL
    rho_sampled_labels <- NULL
  } else {
    if (is_verbose) {
      cat("performing rho sampling...\n")
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
  }
  extinct_labels <- tip_labels[!extant_mask]
  num_extinct <- length(extinct_labels)
  ## We can select which extinctions are deaths, samples and occurrences by
  ## binomially sampling with the correct probabilities.
  num_observed <- rbinom(
    n = 1,
    size = num_extinct,
    prob = params$prob_observed
  )
  num_sampled <- rbinom(
    n = 1,
    size = num_observed,
    prob = params$prob_sampled_given_observed
  )

  observed_labels <- sample(
    x = extinct_labels,
    size = num_observed,
    replace = FALSE
  )
  sampling_labels <- sample(
    x = observed_labels,
    size = num_sampled,
    replace = FALSE
  )

  occurrence_labels <- setdiff(observed_labels, sampling_labels)
  num_occurrences <- length(occurrence_labels)
  unobserved_labels <- setdiff(extinct_labels, observed_labels)

  ## It is useful to have a vector describing what happened to each of the tips
  ## so we will generate this.
  outcome <- character(num_tips)
  for (ix in tip_ix) {
    tl <- tip_labels[ix]
    outcome[ix] <-
      if (is.element(tl, extant_labels)) {
        if (is.null(rho_sampled_labels)) {
          "extant"
        } else if (is.element(tl, rho_sampled_labels)) {
          "rho"
        } else {
          "extant"
        }
      } else if (is.element(tl, unobserved_labels)) {
        "death"
      } else if  (is.element(tl, sampling_labels)) {
        "sampling"
      } else if (is.element(tl, occurrence_labels)) {
        "occurrence"
      } else {
        stop("unaccounted for tip: ", tl)
      }
  }
  ## We then drop the lineages that did not result in a sampled lineage to get
  ## the reconstructed tree.
  rt_tip_labels <- c(sampling_labels, rho_sampled_labels)
  num_rt_tips <- length(rt_tip_labels)
  if (num_rt_tips < 2) {
    stop(paste0(c("\n\n", rep("-", 60), "\n\n\tSIMULATION HAS LESS THAN TWO SEQUENCED SAMPLES!\n\n", rep("-", 60), "\n"), collapse=""))
  }
  rt <- ape::keep.tip(phy, rt_tip_labels)
  rt_tip_and_node_depths <- ape::node.depth.edgelength(rt)
  rt_tip_depths <- head(rt_tip_and_node_depths, num_rt_tips)
  rt_node_depths <- tail(rt_tip_and_node_depths, -num_rt_tips)
  ## The depths of the nodes and tips is relative to the TMRCA so we need to
  ## adjust for this when computing the extact times that they occurred at.
  true_first_sample_time <-
    min(tip_times[outcome == "sampling" | outcome == "rho"])
  depth_first_in_rt <- min(rt_tip_depths)
  rt_offset <- true_first_sample_time - depth_first_in_rt
  ## we can put all of this information into a single dataframe which is more
  ## convenient for subsequent usage.
  event_times_df <- data.frame(
    time = c(-tmrca,
             tip_times[outcome == "occurrence"],
             tip_times[outcome == "sampling"],
             tip_times[outcome == "rho"],
             rt_node_depths + rt_offset
             ),
    event = c("origin",
              rep(
                c("occurrence", "sampling", "rho", "birth"),
                times = c(num_occurrences,
                          num_sampled,
                          ifelse(is.null(num_rho_sampled), 0, num_rho_sampled),
                          ifelse(is.null(num_rho_sampled), num_sampled - 1, num_rt_tips - 1))
              )
              )
  )
  ## Compute the prevalence of infection at the end of the simulation. If there
  ## is no rho sampling this is just the number of extant lineages otherwise we
  ## need to subtract the number that were rho sampled from this.
  if (is.null(params$rho)) {
    final_prevalence <- num_extant
  } else {
    final_prevalence <- num_extant - num_rho_sampled
  }
  if (is_verbose) {
    cat("checking output from run_simulation...\n")
  }
  ## We can do a quick check of the results to make sure that some invariants
  ## have been preserved. Note that one of the checks here will fail a small
  ## percent of the time because it is not exact.
  dur_1 <- max(tip_times[outcome == "sampling" | outcome == "rho"]) + tmrca
  dur_2 <- max(rt_tip_depths) + rt_offset + tmrca
  ## these should be equal up to numerical error.
  stopifnot(abs(dur_1 - dur_2) < 1e-10 * params$duration)
  dur_3 <- max(tip_times[outcome == "extant" | outcome == "rho"]) + tmrca
  stopifnot(abs(dur_3 - params$duration) < time_eps)
  ## If there is aggregation that needs to be carried out then this needs to be
  ## done now just before the simulation results are returned. This operation is
  ## only defined when there is no rho sampling at the end of the simulation.
  if (is.null(params$rho) && (!is.null(params$occ_agg_times) | !is.null(params$seq_agg_times))) {
    ## There is a little extra book keeping involved because the times are
    ## relative to the TMRCA rather than the origin. To keep things consistent
    ## with the TMRCA relative times we need to adjust the aggregation times.
    ## This is why we need to subtract the TMRCA from the given aggregation
    ## times to get the TMRCA relative aggregation times. We add a column for
    ## size to count the number of individuals that were removed in the mock
    ## scheduled event.
    origin_event_row <- filter(event_times_df, event == "origin")
    origin_event_row$size <- NA
    birth_rows <- filter(event_times_df, event == "birth")
    birth_rows$size <- NA
    if (!is.null(params$seq_agg_times)) {
      sampling_times <- event_times_df |> filter(event == "sampling") |> select(time)
      tmrca_rel_seq_agg_times <- params$seq_agg_times - tmrca
      num_agg_obs <- length(tmrca_rel_seq_agg_times) - 1
      if (max(params$seq_agg_times) < max(sampling_times)) {
        stop("There are sequenced samples after the last given aggregation time. It is unclear how to account for the births related to these sequences so you probably want to include another aggregation point to capture them.")
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
      seq_rows <- event_times_df[event_times_df$event == "sampling", ]
      seq_rows$size <- NA
    }
    if (!is.null(params$occ_agg_times)) {
      occurrence_times <- event_times_df |> filter(event == "occurrence") |> select(time)
      tmrca_rel_occ_agg_times <- params$occ_agg_times - tmrca
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
    num_observed = num_observed,
    num_occurrences = num_occurrences))
}

write_plot <- function(simulation_results,
                       parameters,
                       output_directory,
                       is_verbose) {
  hist_plt_df <- data.frame(
    outcome = c("death",
                "sampling",
                "occurrence"),
    empirical =
      c(simulation_results$num_extinct - simulation_results$num_observed,
        simulation_results$num_sampled,
        simulation_results$num_occurrences),
    theory =
      c(qbinom(p = 0.5,
               size = simulation_results$num_extinct,
               prob = 1 - parameters$sampling_prob - parameters$occurrence_prob),
        qbinom(p = 0.5,
               size = simulation_results$num_extinct,
               prob = parameters$sampling_prob),
        qbinom(p = 0.5,
               size = simulation_results$num_extinct,
               prob = parameters$occurrence_prob)),
    theory_min =
      c(qbinom(p = 0.025,
               size = simulation_results$num_extinct,
               prob = 1 - parameters$sampling_prob - parameters$occurrence_prob),
        qbinom(p = 0.025,
               size = simulation_results$num_extinct,
               prob = parameters$sampling_prob),
        qbinom(p = 0.025,
               size = simulation_results$num_extinct,
               prob = parameters$occurrence_prob)),
    theory_max =
      c(qbinom(p = 0.975,
               size = simulation_results$num_extinct,
               prob = 1 - parameters$sampling_prob - parameters$occurrence_prob),
        qbinom(p = 0.975,
               size = simulation_results$num_extinct,
               prob = parameters$sampling_prob),
        qbinom(p = 0.975,
               size = simulation_results$num_extinct,
               prob = parameters$occurrence_prob))
  )

  if (!is.null(parameters$rho)) {
    rho_row <- data.frame(
      outcome = "rho",
      empirical = simulation_results$num_rho_sampled,
      theory = qbinom(p = 0.5,
                      size = simulation_results$num_extant,
                      prob = parameters$rho),
      theory_min = qbinom(p = 0.025,
                          size = simulation_results$num_extant,
                          prob = parameters$rho),
      theory_max = qbinom(p = 0.975,
                          size = simulation_results$num_extant,
                          prob = parameters$rho)
    )
    hist_plt_df <- rbind(hist_plt_df, rho_row)
  }
  tmp <- data.frame(outcome = "extant",
                    empirical = simulation_results$num_extant,
                    theory = NA,
                    theory_min = NA,
                    theory_max = NA)
  hist_plt_df <- rbind(hist_plt_df, tmp)
  hist_plt_df$outcome <- factor(hist_plt_df$outcome, levels = c("death", "occurrence", "sampling", "extant"))
  plt5 <- ggplot(hist_plt_df) +
    geom_col(mapping = aes(x = outcome, y = empirical)) +
    geom_point(mapping = aes(x = outcome, y = theory)) +
    geom_errorbar(mapping = aes(x = outcome, ymin = theory_min, ymax = theory_max)) +
    labs(y = "Count", x = NULL) +
    theme_classic()

  plt6 <- ggplot(data = hist_plt_df) +
    geom_col(mapping = aes(x = outcome, y = empirical), colour = "#636363", fill = "#bdbdbd") +
    scale_x_discrete(NULL,
                     labels = c(
                       "extant" = "Prevalence",
                       "death" = "Unobserved",
                       "sampling" = "Sequenced",
                       "occurrence" = "Unsequenced")) +
    labs(y = NULL, x = NULL) +
    theme_classic()

  is_observed <- is.element(simulation_results$outcome, c("sampling", "occurrence", "rho"))

  tip_annotations <- tibble(
    node = simulation_results$tip_ix,
    outcome = simulation_results$outcome,
    is_observed = is_observed
  )

  tr <- tidytree::treedata(phylo = simulation_results$phylo, data = tip_annotations)

  plt4 <- ggplot(tr, mapping = aes(x, y)) +
    geom_tippoint(mapping = aes(colour = outcome, shape = is_observed),
                  size = 3) +
    geom_nodepoint() +
    geom_vline(data = simulation_results$event_times_df, aes(xintercept = time, colour = event)) +
    geom_tree() +
    labs(colour = "Tip outcome",
         shape = "Observed") +
    theme_tree2(legend.position = "top")

  plt7 <- ggplot(tr, mapping = aes(x, y)) +
    geom_tree(alpha = 0.3) +
    geom_tippoint(mapping = aes(colour = is_observed),
                  size = 2.5) +
    scale_colour_manual(values = c("#bdbdbd", green_hex_colour), labels = c("TRUE" = "Observed", "FALSE" = "Not observed")) +
    labs(colour = NULL) +
    theme_tree(legend.position = "top")

  #' TODO This figure needs fixing because it does not represent rho-samples
  #' properly in the histogram!
  if (is_verbose) {
    cat("writing visualistion to file...\n")
  }
  ggsave(
    filename = paste(
      output_directory,
      "ape-simulation-figure.png",
      sep = "/"
    ),
    plot = plot_grid(plt4, plt5, nrow = 1, rel_widths = c(2, 1)),
    width = 25,
    height = 15,
    units = "cm"
  )
  saveRDS(object = list(plt6, plt7),
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

#' Parse a string \"FROM TO BY\" into a linear vector of values.
parse_from_to_by <- function(from_to_by_string) {
  tmp <- as.numeric(unlist(strsplit(x = from_to_by_string, split = " ")))
  stopifnot(length(tmp) == 3)
  return(seq(from = tmp[1], to = tmp[2], by = tmp[3]))
}

#' Predicate for the command line arguments being valid.
arguments_are_valid <- function(args) {
  if (args$rho > -1.0 && args$seq_agg_times != "") {
    return(FALSE)
  }
  ## We have some assumptions about the input parameters so we should check that
  ## they are satisfied.
  tmp_params <- jsonlite::read_json(args$parameters)
  expected_params <- c("birthRate",
                       "deathRate",
                       "samplingRate",
                       "occurrenceRate")
  for (p in expected_params) {
    if (!is.element(p, names(tmp_params))) {
      warning("Missing parameter: ", p)
    }
  }
  ## If we are being asked to simulate sequences, then there should be a
  ## substitution rate specified in the parameters JSON.
  if (args$simulate_sequences) {
    if (!is.element("substitutionRate", names(tmp_params))) {
      warning("--simulate-sequences flag was given but there is no substitutionRate in the parameters JSON.")
      return(FALSE)
    }
  }
  return(TRUE)
}

main <- function(args) {
  if (args$version) {
    cat(paste(c("ape-sim-",
                paste(as.character(VERSION), collapse="."),
                "\n"), collapse=""))
    return(0)
  }

  if (!arguments_are_valid(args)) {
    print(args)
    stop("The command line arguments are not valid.")
  }

  if (args$verbose) {
    cat("reading parameters from", args$parameters, "\n")
  }
  set.seed(args$seed)

  ## Awkward because you cannot assign NULL via ifelse.
  if (args$rho == -1.0) {
    maybe_rho <-  NULL
  } else {
    maybe_rho <-  args$rho
  }

  ## If there is a request to aggregate the values then we need to parse the
  ## specification of the times at which the aggregation occurs.
  if (args$seq_agg_times == "") {
    maybe_seq_from_to_by <- NULL
  } else {
    maybe_seq_from_to_by <- args$seq_agg_times
  }

  if (args$occ_agg_times == "") {
    maybe_occ_from_to_by <- NULL
  } else {
    maybe_occ_from_to_by <- args$occ_agg_times
  }

  if (args$verbose) {
    if (is.null(maybe_rho)) {
      cat("without rho sampling at the end of the simulation\n")
    } else {
      cat("with rho sampling at the end of the simulation\n")
    }
  }
  params <- read_parameters(
    args$parameters,
    args$duration,
    maybe_rho,
    maybe_seq_from_to_by,
    maybe_occ_from_to_by,
    args$simulate_sequences,
    args$verbose)

  ## the default simulator can produce simulations where the tree is not valid
  ## (because it only has a single tip) so we use the conditioned simulator to
  ## repeat the simulation if this happens.
  sim_result <- run_conditioned_simulation(params, args$verbose)

  if (file.access(args$output_directory, mode = 2) != 0) {
    stop("Cannot write to output directory: ", args$output_directory)
  } else {
    output_filepath <- function(n) {
      paste(args$output_directory, n, sep = "/")
    }
    ## record the events that happened in the simulation
    if (args$verbose) {
      cat("writing output to csv...\n")
    }
    write.table(x = sim_result$event_times_df,
                file = output_filepath("ape-sim-event-times.csv"),
                sep = ",",
                row.names = FALSE)

    if (args$write_newick) {
      if (args$verbose) {
        cat("writing the newick for the reconstructed tree...\n")
      }
      write.tree(phy = sim_result$reconstructed_tree,
                 file = output_filepath("ape-sim-reconstructed-tree.newick"))
    }
    if (args$simulate_sequences) {
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

  if (args$make_plots) {
    write_plot(sim_result, params, args$output_directory, args$verbose)
    if (!is.null(sim_result$aggregated_event_times_df)) {
      write_aggregated_plot(sim_result, params, args$output_directory, args$verbose)
    }
  }
}

if (!interactive()) {
  args <- parser$parse_args()
  main(args)
} else {
  args <- list(
    verbose = TRUE,
    seed = 2,
    parameters = "my-params.json",
    duration = 0.05,
    output_directory = "out",
    make_plots = TRUE,
    rho = -1.0,
    ## rho = 0.5,
    seq_agg_times = "",
    occ_agg_times = "",
    write_newick = FALSE,
    simulate_sequences = TRUE
  )
  main(args)
}
