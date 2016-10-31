#' Generate a random SIR simulation
#'
#' @return
#' @export
#'
#' @examples
example_sim_sir <- function () {
  simulate_sir(c(R0=2, Tg=5, k=1, N=5000, S=4999, I=1), dt=0.1,
               max_num_dt=1000, min_epi_size=20, maximum_attempts=1,
               track_transmissions=TRUE)
}
