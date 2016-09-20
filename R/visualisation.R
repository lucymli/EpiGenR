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
