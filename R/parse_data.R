count_cases <- function (times, time_vector) {
  sapply(time_vector, function (time_step) {
    sum(times == time_step)
  })
}

#' Generate a time-series based on line list data
#'
#' @param line_list A matrix or dataframe in which the first column is the ID number of each case and the second column is the time of infection.
#' @param step_size Size of time step. Cases occurring within each time step are accumulated
#' @param split_by A factor that splits the line_list. The output would thus be length(levels(split_by)) number of time-series.
#'
#' @return
#' @export
#'
#' @examples
time_series_from_line_list <- function (line_list, step_size=1, keep.dates=FALSE, split_by=NULL) {
  if (is.list(line_list)) {
    if ("total_sampled" %in% names(line_list)) line_list <- line_list$infected_sampled[, 3]
    else line_list <- line_list$infected[, 3]
  }
  if (keep.dates) {
    line_list_dates <- line_list
    min_date <- min(line_list_dates)
    if (!is.null(split_by)) {
      line_list_dates <- split(line_list, split_by)
      min_date <- lapply(line_list_dates, min)
    }
  }
  if (class(line_list)=="Date") {
    line_list <- as.numeric(line_list - min(line_list))
  }
  rounded_times <- ceiling(line_list/step_size)
  max_rounded_time <- max(rounded_times)
  time_vec <- 1:max_rounded_time
  if (is.null(split_by)) {
    counts <- count_cases(rounded_times, time_vec)
    time_series <- data.frame(time=time_vec*step_size, incidence=counts)
  } else {
    incidence <- lapply(split(rounded_times, split_by), function (x) {
      counts <- count_cases(x, time_vec)
      return(counts)
    })
    time_series <- data.frame(time=time_vec*step_size, do.call(cbind, incidence))
    names(time_series)[-1] <- as.character(split_by)
  }
  if (keep.dates) {
    if (is.null(split_by)) {
      time_series[, 1] <- min_date + (seq_len(nrow(time_series))-1)*step_size
    } else {
      time_series <- lapply(time_series, function (x) {
        x[, 1] <- min_date + (seq_len(nrow(x))-1)*step_size
        x
      })
    }
  }
  return(time_series)
}



#' Sample from a time series
#'
#' @param time_series A matrix with 2 columns: time and count
#' @param prob A numeric denoting the probability of sampling
#' @param type A string indicating the sampling strategy.
#'
#' @details There are 3 sampling strategies: 1. uniform, 2. proportional and 3. log-proportional.
#' In uniform sampling, a fixed number of samples is taken at each time step. In proportional
#' sampling, the number of samples taken is proportional to the size of the population at
#' each time step. In log-proportional sampling, the number of samples is proportional to
#' the log of the population size at each time step, which could happen during an outbreak
#' if a laboratory is overwhelmed with samples. In all cases, the ratio of the resulting
#' number of samples and the total population size should be roughly equal to prob.
#'
#' @return
#' @export
#'
#' @examples
sampled_from_time_series <- function (time_series, prob, type="uniform") {
  total <- sum(time_series[, 2])
  sampled_time_series <- cbind(time_series[, 1], 0)
  time_steps_with_counts <- which(time_series[, 2] > 0)
  total_sampled <- rbinom(1, total, prob)
  if (type%in%c("uniform", "log-proportional")) {
    if (type == "uniform") PROB <- NULL
    else PROB <- log(time_series[time_steps_with_counts, 2])
    sampled_time_steps <- table(sample(time_steps_with_counts, total_sampled, prob=PROB, replace=TRUE))
    for (i in seq_along(sampled_time_steps)) {
      time_step <- as.numeric(names(sampled_time_steps)[i])
      sampled_time_series[time_step, 2] <- min(sampled_time_steps[i], time_series[time_step, 2])
    }
  }
  if (type == "proportional") {
    sampled_time_series[, 2] <- rbinom(nrow(time_series), time_series[, 2], prob)
  }
  return(sampled_time_series)
}


get_tip_sample_times <- function (tr) {
  tip_labels <- tr$tip.label
  tip_labels_mat <- do.call(rbind, strsplit(tip_labels, "_"))
  if (any(grepl("-", tip_labels_mat[, 2]))) {
    tip_sample_times <- lubridate::decimal_date(as.Date(tip_labels_mat[, 2]))
  } else {
    tip_sample_times <- as.numeric(tip_labels_mat[, 2])
  }
  return(tip_sample_times)
}

#' Obtain phylogenetic data in the form of a time-series
#'
#' @param tr
#' @param step_size
#'
#' @return
#' @export
#'
#' @examples
time_series_from_tree <- function (tr, step_size=1) {
  time_series_from_line_list(get_tip_sample_times(tr), step_size=step_size)
}

get_last_tip_time <- function (tr) {
  return(max(get_tip_sample_times(tr)))
}

#' Aligns the epidemic time-series with the genealogical time-series in time
#'
#' @param epi A data frame where the first column denotes time and the second column onwards contains the count
#' @param gen Genealogical time series
#' @param dt Size of time step
#' @param last_tip_time A numeric value denoting the time that the most recent tip was sampled
#'
#' @return
#' @export
#'
#' @examples
align_epi_gen_data <- function (epi, gen, dt, last_tip_time) {
  num_dt_gen <- length(gen)
  last_tip_dt <- ceiling(last_tip_time/dt)
  tmrca_dt <- last_tip_dt - num_dt_gen + 1
  first_epi_dt <- round(epi[1, 1]/dt)
  last_epi_dt <- round(epi[nrow(epi), 1]/dt)
  if (tmrca_dt > first_epi_dt) {
    gen <- c(lapply(1:(tmrca_dt-first_epi_dt), function (x) list(binomial=0, intervals=0)), gen)
  }
  if (last_tip_dt < last_epi_dt) {
    gen <- c(gen, lapply(1:(last_epi_dt-last_tip_dt), function (x) list(binomial=0, intervals=0)))
  }
  if (tmrca_dt < first_epi_dt) {
    epi <- rbind(cbind(seq(tmrca_dt*dt, by=dt, length.out=first_epi_dt-tmrca_dt),
                       matrix(0, nrow=first_epi_dt-tmrca_dt, ncol=ncol(epi)-1)), epi)
  }
  if (last_tip_dt > last_epi_dt) {
    epi1 <- cbind(seq(to=last_tip_dt*dt, by=dt, length.out=last_tip_dt- last_epi_dt),
                            matrix(0, nrow=last_tip_dt- last_epi_dt, ncol=ncol(epi)-1))
    if (is.data.frame(epi)) {
      epi1 <- data.frame(epi1)
      names(epi1) <- names(epi)
    }
    epi <- rbind(epi, epi1)
  }
  return(list(epi=epi, gen=gen))
}

#' Produce data from line list, phylogeny, or both.
#'
#' @param epi
#' @param phy
#' @param dt
#'
#' @return
#' @export
#'
#' @examples
get_data <- function (epi=NULL, phy=NULL, dt=1) {
  epi_data <- gen_data <- NULL
  use.epi <- !is.null(epi)
  use.gen <- !is.null(phy)
  if (!is.null(epi)) {
    epi_data <- time_series_from_line_list(epi, dt)
  }
  if (!is.null(phy)) {
    gen_data <- time_series_from_tree(phy, dt)
  }
  if (use.epi&&use.gen) {
    all_data <- align_epi_gen_data(epi_data, gen_data, dt, get_last_tip_time(phy))
  } else if (use.epi) {
    all_data <- epi_data
  } else {
    all_data <- gen_data
  }
  return (all_data)
}
