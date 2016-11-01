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


#' Plot the posterior distribution of pMCMC analysis with one, or both of epi and phylogenetic data
#'
#' @param param.name
#' @param dataset
#' @param new.param.name
#' @param fill.colours
#' @param show.legend
#'
#' @return
#' @export
#'
#' @examples
plot.posterior <- function (param.name, dataset, new.param.name=NULL, fill.colours=NULL, show.legend=TRUE) {
  require (ggplot2)
  require(lubridate)
  require(coda)
  hpd.interval <- do.call(data.frame, lapply(split(dataset[, param.name], dataset$data.type), function (x) {
    if (is.Date(x[1])) {
      parsed.x <- as.numeric(x-(min(x)-1))
      hpd.interval <- hpd(parsed.x) + (min(x)-1)
      #hpd.interval <- as.Date(date_decimal(EpiGenR::hpd(decimal_date(unlist(x)))))
    }
    else hpd.interval <- EpiGenR::hpd(x)
    return (hpd.interval)
  }))
  Plot <- ggplot(dataset) +
    theme_bw() +
    geom_density(aes_string(x=param.name, fill="data.type"), alpha=.4) +
    scale_fill_manual(name="Data")
  if (is.Date(hpd.interval[1,1])) {
    Plot <- Plot +
      scale_x_date(limits=c(min(do.call("c", hpd.interval[2, ])),
                            max(do.call("c", hpd.interval[3, ]))))
  } else {
    Plot <- Plot +
      xlim(min(unlist(hpd.interval[2, ])), max(unlist(hpd.interval[3, ])))
  }
  if (!is.null(fill.colours)) {
    Plot <- Plot + scale_fill_manual(values=fill.colours, name="Data")
  }
  if (!is.null(new.param.name)) {
    Plot <- Plot + xlab(new.param.name)
  }
  if (!show.legend) {
    Plot <- Plot + theme(legend.position="none")
  }
  return (Plot)
}
