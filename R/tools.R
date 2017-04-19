#' Sum a vector of numbers in non-overlapping groups of x
#'
#' @param vec
#' @param every
#'
#' @return
#' @export
#'
#' @examples
sum_every <- function (vec, every = 1) {
  vec <- as.numeric(vec)
  len <- length(vec)
  begin <- seq(1, len, every)
  end <- seq(every, len, every)
  if (length(begin) > length(end))
    end <- c(end, len)
  out <- sapply(1:length(begin), function(i) {
    sum(vec[begin[i]:end[i]])
  })
  return(out)
}


#' Calculates the corresponding lognormal parameters for a given mean and sd
#'
#' @param MEAN
#' @param SD
#'
#' @return
#' @export
#'
#' @examples
get_lognormal_params <- function (MEAN, SD) {
  sigma <- sqrt(log(1 + (SD/MEAN)^2))
  zeta <- log(MEAN/sqrt(1+(SD/MEAN)^2))
  return(c(zeta, sigma))
}


#' Generates random numbers according to the log-normal distribution with mean and sd
#' defined in non transformed space
#'
#' @param N
#' @param MEAN
#' @param SD
#'
#' @return
#' @export
#'
#' @examples
rlnorm_real <- function (N, MEAN, SD) {
  transformed <- get_lognormal_params(MEAN, SD)
  rlnorm(N, transformed[1], transformed[2])
}

#' Get the median and 95% HPD interval of a vector
#'
#' @param vec
#' @param conf
#'
#' @return
#' @export
#'
#' @examples
hpd <- function (vec, conf=0.95, show.median=TRUE) {
  x <- as.numeric(vec)
  if (all(is.na(x))) return (rep(NA, 3))
  med <- median(x, na.rm=TRUE)
  ci <- coda::HPDinterval(mcmc(x), conf)
  out <- c(median=med, lower=ci[1], upper=ci[2])
  if (show.median) return (out)
  else return (out[-1])
}

hpd.Date <- function (vec, conf=0.95, show.median=TRUE) {
  if (all(is.na(vec))) return (rep(NA, 3))
  start.date <- min(vec[!is.na[vec]])-1
  x <- as.numeric(vec[!is.na[vec]] - start.date)
  med <- median(x, na.rm=TRUE)
  ci <- coda::HPDinterval(coda::mcmc(x), conf)
  out <- c(median=med, lower=ci[1], upper=ci[2]) + start.date
  if (show.median) return (out)
  else return (out[-1])
}

str_to_sentence <- function (strings) {
  if (length(strings)==1) return (strings)
  if (length(strings)==2) return (paste(strings, collapse=" and "))
  part1 <- paste(strings[1:(length(strings)-1)], collapse=", ")
  return(paste(part1, tail(strings, 1), sep=", and "))
}

values_ci_to_str <- function (values, sep=", ") {
  paste0(values[1], " (", values[2], sep, values[3], ")")
}
