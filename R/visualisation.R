#' Extracts the legend from a plot
#'
#' @param A ggplot object
#'
#' @return
#' @export
#'
#' @examples
extract_gglegend <- function (Plot) {
  legend.tmp <- ggplot_gtable(ggplot_build(Plot))
  legend.plot <- legend.tmp$grobs[[which("guide-box" == sapply(legend.tmp$grobs, `[[`, "name"))]]
  return (legend.plot)
}

#' Plot the MCMC trace using the results data frame
#'
#' @param x data frame where each row contains the parameter values at a particular MCMC iteration
#' @param burnin If not NULL, the proportion of the chain to discard.
#' @param excludes A vector of character. Column/variable names in the data frame x to exclude from visualisation
#'
#' @return
#' @export
#'
#' @examples
plot.mcmc.trace <- function (x, burnin=NULL, excludes=NULL) {
  if (!is.null(burnin)) {
    x <- x[round(nrow(x)*burnin):nrow(x), ]
  }
  if (!is.null(excludes)) {
    y <- x[, !(names(x) %in% excludes)]
  }
  y$iteration <- 1:nrow(y)
  y <- reshape2::melt(y, id.vars="iteration")
  require(ggplot2)
  ggplot(y) + geom_line(aes(x=iteration, y=value)) +
    facet_grid(variable~., scales="free_y")
}
