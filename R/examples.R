#' Generate a random SIR simulation
#'
#' @return
#' @export
#'
#' @examples
example.sim.sir <- function () {
  simulate_sir(c(R0=2, Tg=5, k=1, N=5000, S=4999, I=1), 0.1, 1000, 20, 1, FALSE)
}
