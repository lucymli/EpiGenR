#' Generate the input files for the C++ implementation of pMCMC
#'
#' @param params  Model parameters
#' @param mcmc_options MCMC options
#' @param data Data
#'
#' @return
#' @export
#'
#' @examples
generate_cpp_input_files <- function(dt=NULL, params=NULL, mcmc_options=NULL, initial_states=NULL, data=NULL,
                                     params_file=NULL, mcmc_options_file=NULL, initial_states_file=NULL, data_file=NULL) {
  if (is.null(params_file)) {
    params_file <- tempfile(fileext="_params.txt")
  }
  if (!is.null(params)) {
    writeLines(unlist(lapply(params, paste, collapse=" ")), params_file)
  }
  if (is.null(mcmc_options_file)) {
    mcmc_options_file <- tempfile(fileext="_mcmc_options.txt")
    file.create(mcmc_options_file)
  }
  if (!is.null(mcmc_options)) {
    out <- sapply(seq_along(mcmc_options), function (i) {
      cat(names(mcmc_options)[i], " ", mcmc_options[i], "\n", sep="", file=mcmc_options_file, append=TRUE)
    })
  }
  if (is.null(initial_states_file)) {
    initial_states_file <- tempfile(fileext="_initial_states.txt")
  }
  if (!is.null(initial_states)) {
    cat(initial_states, sep="\n", file=initial_states_file)
  }
  if (is.null(data_file)) {
    data_file <- tempfile()
    if ((length(data)==2) && is.list(data)) which_lik <- 0
    else if (is.data.frame(data) || is.matrix(data)) which_lik <- 1
    else which_lik <- 2
    if (which_lik==0) num_dt <- nrow(data[[1]])
    else if (which_lik==1) num_dt <- nrow(data)
    else num_dt <- length(data)
    filenames <- c()
    if (which_lik<2) {
      epi_name <- paste0(data_file, "_epi_data.txt")
      if (which_lik==0) epi_data <- data[[1]]
      else epi_data <- data
      if (ncol(epi_data) == 2) {
        cat(num_dt, dt, 1, epi_data[, 2], sep="\n", file=epi_name)
      } else {
        cat(num_dt, dt, ncol(epi_data)-1, sep="\n", file=epi_name)
        write.table(epi_data[, -1], quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
                    file=epi_name, append=TRUE)
      }
      filenames <- c(filenames, epi_name)
    }
    if (which_lik!=1) {
      mcmc_options["num_trees"] <- 1
      gen_name <- paste0(data_file, "_gen_data.txt")
      cat(num_dt, dt, sep="\n", file=gen_name)
      if (which_lik==0) gen_data <- data[[2]]
      else gen_data <- data
      out <- lapply(gen_data, function (x) {
        cat(paste(x$binomial, collapse=" "), sep="\n", file=gen_name, append=TRUE)
      })
      out <- lapply(gen_data, function (x) {
        cat(paste(x$intervals, collapse=" "), sep="\n", file=gen_name, append=TRUE)
      })
      filenames <- c(filenames, gen_name)
    }
    data_file <- filenames
  }
  cpp_input_files <- c(params_file=params_file, mcmc_options_file=mcmc_options_file,
                       initial_states_file=initial_states_file, data_files=paste(data_file, collapse=" "))
  commandline.command <- paste(cpp_input_files, collapse=" ")
  return (commandline.command)
}

create_mcmc_options <- function (particles=1000, iterations=1000, log_every=1, pfilter_every=20,
                                 which_likelihood=0, pfilter_threshold=1.0, num_threads=4, heat_factor=NA,
                                 heat_length=NA, cool_rate=NA, log_filename="log.txt", traj_filename="traj.txt",
                                 save_traj=NA,
                                 use_lhs=NA, lhs_divides=NA, lhs_iterations=NA) {
  mcmc_options <- c(particles=particles, iterations=iterations, log_every=log_every,
                    pfilter_every=pfilter_every, which_likelihood=which_likelihood,
                    pfilter_threshold=pfilter_threshold, num_threads=num_threads,
                    heat_factor=heat_factor, heat_length=heat_length, cool_rate=cool_rate,
                    log_filename=log_filename,
                    traj_filename=traj_filename, save_traj=save_traj,
                    use_lhs=use_lhs, lhs_divides=lhs_divides, lhs_iterations=lhs_iterations)
  return (mcmc_options[!is.na(mcmc_options)])
}


create_params_list <- function (param_names=c("param"), init_param_values=c(1),
                                params_to_estimate=c("param"), transform=c(NA),
                                prior=c("unif"), prior_params=list(c(0.0, 0.1)),
                                proposal=NULL, proposal_params=list(c(0.1, 0.0, 1.0)),
                                optimal_acceptance=0.234, lower_acceptance=0.1, upper_acceptance=0.8,
                                adapt_every=20, max_adapt_times=100) {
  estimate <- param_names %in% params_to_estimate
  num_params <- length(param_names)
  input_params1 <- list(c(total_params=num_params, optimal_acceptance=optimal_acceptance,
                         lower_acceptance=lower_acceptance, upper_acceptance=upper_acceptance,
                         adapt_every=adapt_every, max_adapt_times=max_adapt_times))
  input_params2 <- data.frame(matrix(nrow=num_params, ncol=12))
  names(input_params2) <- c("param_value", "param_name", "estimate", "transform",
                            "prior", "prior_param1", "prior_param2", "prior_param3",
                            "proposal", "proposal_sd", "proposal_min", "proposal_max"
                            )
  input_params2[, c(paste0("prior_param", 1:3), "proposal_sd", "proposal_min", "proposal_max")] <- 0.0
  input_params2$param_value <- init_param_values
  input_params2$param_name <- param_names
  input_params2$estimate <- tolower(estimate)
  input_params2$transform <- NA
  input_params2$prior <- NA
  input_params2$prior[estimate] <- prior
  input_params2$proposal <- NA
  if (is.null(proposal)) input_params2$proposal[estimate] <- "normal"
  a <- lapply(1:length(prior_params), function (i) {
    j <- which(estimate)[i]
    input_params2$transform[j] <<- transform[i]
    input_params2[j, "prior_param1"] <<- prior_params[[i]][1]
    if (length(prior_params[[i]]) > 1) input_params2[j, "prior_param2"] <<- prior_params[[i]][2]
    if (length(prior_params[[i]]) > 2) input_params2[j, "prior_param3"] <<- prior_params[[i]][3]
    input_params2$proposal_sd[j] <<- proposal_params[[i]][1]
    input_params2$proposal_min[j] <<- proposal_params[[i]][2]
    input_params2$proposal_max[j] <<- proposal_params[[i]][3]
  })
  params_list <- c(input_params1, lapply(1:nrow(input_params2), function (i) input_params2[i, ]))
  return(params_list)
}


#' Run the compile C++ program with the given input
#'
#' @param program
#' @param input
#' @param wait
#'
#' @return
#' @export
#'
#' @examples
run_pMCMC <- function (program, input, wait) {
  out <- system(paste(program, input), intern=wait, wait=wait)
  return (out)
}
