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
