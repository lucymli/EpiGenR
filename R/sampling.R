downsample <- function (epidemic_obj, strategy="proportional", prob=NULL, num.samples=NULL, time.step=NULL) {
  if (strategy=="proportional") {
    if (is.null(prob)) prob <- 1
    selected <- which(rbinom(nrow(epidemic_obj$infected), 1, prob)==1)
  } else {
    if (is.null(num.samples)) num.samples <- 1
    if (is.null(time.step)) time.step <- 1
    if (strategy=="log-proportional") {
      time.splits <- seq(0, max(epidemic_obj$infected[, 3]), length.out=min(nrow(epidemic_obj$total_dt), 100))
    } else {
      time.splits <- seq(0, max(epidemic_obj$infected[, 3]), by=time.step)
    }
    selected <- unlist(lapply(1:(length(time.splits)-1), function (i) {
      selected <- which(findInterval(epidemic_obj$infected[, 3], time.splits[c(i, i+1)]) == 1)
      if (length(selected)==0) return (NULL)
      if (strategy=="log-proportional") {
        sample.size <- rbinom(1, log(length(selected)), prob)
      } else {
        sample.size <- min(length(selected), num.samples)
      }
      indices <- sample(selected, sample.size)
      return (indices)
    }))
  }
  selected <- selected[epidemic_obj$infected[selected, 1]>=0]
  subepidemic <- epidemic_obj$infected[selected, ]
  subtotal <- nrow(subepidemic)
  epidemic_obj$infected_sampled <- subepidemic
  epidemic_obj$total_sampled <- subtotal
  epidemic_obj$sampled_individuals <- selected
  return (epidemic_obj)
}
