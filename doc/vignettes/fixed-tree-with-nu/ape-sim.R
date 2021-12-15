#!/usr/bin/env Rscript
#'
VERSION <- c(0,1,4)
#' ape-sim-0.1.4
#' =============
#'
#' Use the ape package to simulate the BDSCOD process from the command line.
#'
#' Simulations are conditioned on there being more than one sequenced sample. If
#' no satisfactory simulation is generated in 100 replicates then the program
#' will crash. The parameters for the simulation are read in from an XML file as
#' demonstrated below.
#'
#' Usage
#' -----
#'
#' The following will simulate without either a catastrophe or a disaster:
#'
#' $ ./ape-sim.R --verbose demo.xml
#'
#' where the \code{demo.xml} file should like the following
#'
#' <?xml version="1.0" encoding="UTF-8" standalone="no"?>
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
suppressPackageStartupMessages(library(XML))

green_hex_colour <- "#1b9e77"
purple_hex_colour <- "#7570b3"

#' Parse the XML configuration of the simulation. See the documentation above
#' for an example.
parse_xml_configuration <- function(filepath) {
  xml_input <- xmlToList(xmlParse(filepath))
  xml_opts <- as.list(xml_input$configuration$options)
  xml_params <- as.list(xml_input$configuration$parameters)

  params = list(
    birth_rate = as.numeric(xml_params$birthRate),
    death_rate = as.numeric(xml_params$deathRate),
    sampling_rate = as.numeric(xml_params$samplingRate),
    occurrence_rate = as.numeric(xml_params$occurrenceRate),
    duration = as.numeric(xml_params$duration)
  )

  if (is.element("rho", names(xml_params))) {
    params$rho <- as.numeric(xml_params$rho)
  } else {
    params$rho <- NULL
  }

  if (is.element("nu", names(xml_params))) {
    params$nu <- as.numeric(xml_params$nu)
  } else {
    params$nu <- NULL
  }

  if (is.element("substitutionRate", names(xml_params))) {
    params$subsitution_rate <- as.numeric(xml_params$substitutionRate)
  }

  if (is.element("seqLength", names(xml_params))) {
    params$seq_length <- as.numeric(xml_params$seqLength)
  }

  params$net_rem_rate <- params$death_rate + params$sampling_rate + params$occurrence_rate
  params$net_per_capita_event_rate <- params$birth_rate + params$net_rem_rate
  params$sampling_prob <- params$sampling_rate / params$net_rem_rate
  params$occurrence_prob <- params$occurrence_rate / params$net_rem_rate
  params$prob_observed <- params$sampling_prob + params$occurrence_prob
  params$prob_sampled_given_observed <-
    params$sampling_prob / params$prob_observed

  opts <- list(
    seed = as.integer(xml_opts$seed),
    output_directory = xml_opts$outputDirectory,
    make_plots = as.logical(xml_opts$makePlots),
    write_newick = as.logical(xml_opts$writeNewick),
    simulate_sequences = as.logical(xml_opts$simulateSequences))

  if (is.element("seq_agg_times", names(xml_opts))) {
    opts$seq_agg_times <- parse_from_to_by(xml_opts$seq_agg_times)
  } else {
    opts$seq_agg_times <- NULL
  }

  if (is.element("occ_agg_times", names(xml_opts))) {
    opts$occ_agg_times <- parse_from_to_by(xml_opts$occ_agg_times)
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
         help = "Filepath to XML configuration")
## parser$add_argument(
##          "--seq-agg-times",
##          type = "character",
##          default = "",
##          help = "Specification of aggregation times for sequenced samples: \"FROM TO BY\". These values get read as three numbers then form the arguments for the seq function.")
## parser$add_argument(
##          "--occ-agg-times",
##          type = "character",
##          default = "",
##          help = "Specification of aggregation times for unsequenced samples (occurrence data). See the details of --seq-agg-times.")
## parser$add_argument(
##          "--simulate-sequences",
##          action = "store_true",
##          default = FALSE,
##          help = "Simulate sequences for each of the leaves of the reconstructed tree. This makes use of the simSeq function from phangorn. If this parameter is given, then it is assumed that there will be a substitutionRate given in the parameters JSON.")

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

sequence_simulation <- function(tr, len, sub_rate) {
  return(simSeq(tr, l = len, rate = sub_rate))
}

run_simulation <- function(params, options, is_verbose) {
  time_eps <- 1e-10
  if (is_verbose) {
    cat("simulating tmrca and phylogeny...\n")
  }
  ## because the trees generated by rlineage start from the TMRCA rather than
  ## the origin we need to simulate the length of the root first and subtract
  ## this from the duration.
  tmrca <- rexp(n = 1, rate = params$net_per_capita_event_rate)
  stopifnot(params$duration > tmrca)
  phy <- rlineage(
    params$birth_rate,
    params$net_rem_rate,
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
  extinct_labels <- tip_labels[!extant_mask]
  num_extinct <- length(extinct_labels)
  ## We need to do a the appropriate sampling at the end of the simulation if
  ## either a rho- or nu-sample has been requested, otherwise we need to
  ## propagate the null values.
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
  ## We can select from the extinct lineages which will be designated as deaths,
  ## samples and occurrences by binomially sampling with the correct
  ## probabilities and then uniformly selecting this many.
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
        if (is.null(c(rho_sampled_labels, nu_sampled_labels))) {
          "extant"
        } else if (is.element(tl, rho_sampled_labels)) {
          "rho"
        } else if (is.element(tl, nu_sampled_labels)) {
          "nu"
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
  ## We then drop the lineages that did not result in a sampled lineage and the
  ## remaining tree is the reconstructed tree.
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
  ## convenient for export and subsequent usage.
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
  ## extant lineages otherwise we need to subtract the number that were
  ## rho-sampled or nu-sampled depending on which was carried out.
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
  ## We can do a quick check of the results to make sure that some invariants
  ## have been preserved. Note that one of the checks here will fail a small
  ## percent of the time because it is not exact.
  dur_0 <- params$duration
  dur_1 <- max(tip_times[outcome == "sampling" | outcome == "rho" | outcome == "nu"]) + tmrca
  dur_2 <- max(rt_tip_depths) + rt_offset + tmrca
  dur_3 <- max(tip_times[outcome == "extant" | outcome == "rho" | outcome == "nu"]) + tmrca
  ## these should be equal up to numerical error.
  fudge <- 1e-3
  if (!(abs(dur_0 - dur_1) < fudge * params$duration)) {
    cat("measure 1 is ", dur_1, "\n")
    stop("measures of simulation duration 1 is too different from requested duration.")
  }
  if (!(abs(dur_0 - dur_2) < fudge * params$duration)) {
    cat("measure 2 is ", dur_2, "\n")
    stop("measures of simulation duration 2 is too different from requested duration.")
  }
  if (!(abs(dur_0 - dur_3) < fudge * params$duration)) {
    cat("measure 3 is ", dur_3, "\n")
    stop("measures of simulation duration 3 is too different from requested duration.")
  }
  if (!(abs(dur_1 - dur_2) < fudge * params$duration)) {
    stop("measures of simulation duration seem too different.")
  }

  ## TODO The aggregation code below should probably be refactored and applied
  ## as a post-simulation step rather than being included here.

  ## If there is aggregation that needs to be carried out then this needs to be
  ## done now just before the simulation results are returned. This operation is
  ## only defined when there are no rho-samples and no nu-samples at the end of
  ## the simulation (which should be enforced by the validation of the
  ## parameters...)
  if (is.null(c(params$rho, params$nu)) &&
      (!is.null(options$occ_agg_times) | !is.null(options$seq_agg_times))) {
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
    if (!is.null(options$seq_agg_times)) {
      sampling_times <- event_times_df |> filter(event == "sampling") |> select(time)
      tmrca_rel_seq_agg_times <- options$seq_agg_times - tmrca
      num_agg_obs <- length(tmrca_rel_seq_agg_times) - 1
      if (max(options$seq_agg_times) < max(sampling_times)) {
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
      ## If there are not any sequence aggregation times then the sequenced
      ## samples still need to be represented by sampling events so we pull
      ## these values out and store them in an appropriate data frame for
      ## subsequent use.
      seq_rows <- event_times_df[event_times_df$event == "sampling", ]
      seq_rows$size <- NA
    }
    if (!is.null(options$occ_agg_times)) {
      occurrence_times <- event_times_df |> filter(event == "occurrence") |> select(time)
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
  if (!is.null(parameters$nu)) {
    nu_row <- data.frame(
      outcome = "nu",
      empirical = simulation_results$num_nu_sampled,
      theory = qbinom(p = 0.5,
                      size = simulation_results$num_extant,
                      prob = parameters$nu),
      theory_min = qbinom(p = 0.025,
                          size = simulation_results$num_extant,
                          prob = parameters$nu),
      theory_max = qbinom(p = 0.975,
                          size = simulation_results$num_extant,
                          prob = parameters$nu)
    )
    hist_plt_df <- rbind(hist_plt_df, nu_row)
  }
  tmp <- data.frame(outcome = "extant",
                    empirical = simulation_results$num_extant,
                    theory = NA,
                    theory_min = NA,
                    theory_max = NA)
  hist_plt_df <- rbind(hist_plt_df, tmp)
  hist_plt_df$outcome <- factor(hist_plt_df$outcome, levels = c("death", "occurrence", "sampling", "extant", "rho", "nu"))
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

  is_observed <- is.element(simulation_results$outcome, c("sampling", "occurrence", "rho", "nu"))

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
  if (from_to_by_string == "") {
    NULL
  } else {
    tmp <- as.numeric(unlist(strsplit(x = from_to_by_string, split = " ")))
    stopifnot(length(tmp) == 3)
    seq(from = tmp[1], to = tmp[2], by = tmp[3])
  }
}

#' Predicate for the configuration read from XML being valid.
configuration_is_valid <- function(config) {
  params <- config$params
  opts <- config$options
  if ((!is.null(params$rho)) && (params$seq_agg_times != "")) {
    return(FALSE)
  }
  ## If we are being asked to simulate sequences, then there should be a
  ## substitution rate specified in the parameters JSON.
  if (opts$simulate_sequences) {
    if (!is.element("substitution_rate", names(params))) {
      warning("--simulate-sequences flag was given but there is no substitutionRate in the parameters XML.")
      return(FALSE)
    }
    if (!is.element("seq_length", names(params))) {
      warning("--simulate-sequences flag was given but there is no seqLength in the parameters XML.")
      return(FALSE)
    }
  }
  print(names(params))
  print(opts)
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
  if (!dir.exists(opts$output_directory)) {
    warning("missing output directory: ", opts$output_directory)
    return(FALSE)
  }
  ## If we are being asked to use a particular output directory it should exist.
  return(TRUE)
}

main <- function(args, config) {
  if (args$version) {
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

  ## This branch demonstrates how to construct an appropriate XML configuration
  ## and simulate from it. This is particularly useful for development but may
  ## be useful as a way to run the code...
  args <- list(
    verbose = TRUE,
    version = FALSE,
    xml = tempfile()
  )

  parameters_node <- xmlNode(
    "parameters",
    attrs = list(birthRate = "3.0",
                 deathRate = "1.0",
                 samplingRate = "0.5",
                 occurrenceRate = "0.5",
                 nu = "0.5",
                 ## rho = "-1.0",
                 duration = "4.0"))

  options_node <- xmlNode(
    "options",
    attrs = list(seed = "1",
                 writeNewick = "true",
                 makePlots = "true",
                 outputDirectory = "out",
                 simulateSequences = "false",
                 seq_agg_times = "",
                 occ_agg_times = ""))

  input_node <- xmlNode(
    "configuration",
    parameters_node,
    options_node
  )

  z <- xmlNode(
    "ape",
    input_node,
    attrs = list(version = "0.1.2")
  )

  saveXML(
    z,
    file = args$xml,
      prefix = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
  )

  config <- parse_xml_configuration(args$xml)
  main(args, config)
}
