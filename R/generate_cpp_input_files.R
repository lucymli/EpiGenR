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
generate_cpp_input_files <- function(dt, params, mcmc_options, initial_states, data,
                                     params_file=NULL, mcmc_options_file=NULL, initial_states_file=NULL, data_file=NULL) {
  if (is.null(params_file)) params_file <- tempfile(fileext="_params.txt")
  if (is.null(mcmc_options_file)) mcmc_options_file <- tempfile(fileext="_mcmc_options.txt")
  if (is.null(initial_states_file)) initial_states_file <- tempfile(fileext="_initial_states.txt")
  if (is.null(data_file)) data_file <- tempfile()
  # file.create(data_file)
  num_dt <- nrow(data[[1]])
  filenames <- c()
  if (mcmc_options["which_likelihood"]<2) {
    epi_name <- paste0(data_file, "_epi_data.txt")
    if (ncol(data[[1]]) == 2) {
      cat(num_dt, dt, 1, data[[1]][, 2], sep="\n", file=epi_name)
    } else {
      cat(num_dt, dt, ncol(data[[1]])-1, sep="\n", file=epi_name)
      write.table(data[[1]][, -1], quote=FALSE, row.names=FALSE, col.names=FALSE, sep=" ",
                  file=epi_name, append=TRUE)
    }
    filenames <- c(filenames, epi_name)
  }
    if (mcmc_options["which_likelihood"]!=1) {
    mcmc_options["num_trees"] <- 1
    gen_name <- paste0(data_file, "_gen_data.txt")
    cat(num_dt, dt, sep="\n", file=gen_name)
    out <- lapply(data[[2]], function (x) {
      cat(paste(x$binomial, collapse=" "), sep="\n", file=gen_name, append=TRUE)
    })
    out <- lapply(data[[2]], function (x) {
      cat(paste(x$intervals, collapse=" "), sep="\n", file=gen_name, append=TRUE)
    })
    filenames <- c(filenames, gen_name)
  }
  writeLines(unlist(lapply(params, paste, collapse=" ")), params_file)
  file.create(mcmc_options_file)
  out <- sapply(seq_along(mcmc_options), function (i) {
    cat(names(mcmc_options)[i], " ", mcmc_options[i], "\n", sep="", file=mcmc_options_file, append=TRUE)
  })
  cat(initial_states, sep="\n", file=initial_states_file)
  cpp_input_files <- c(params_file=params_file, mcmc_options_file=mcmc_options_file, initial_states_file=initial_states_file, data_files=paste(filenames, collapse=" "))
  commandline.command <- paste(cpp_input_files, collapse=" ")
  return (commandline.command)
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
