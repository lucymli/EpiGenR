# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Simulate an epidemic according to an SIR model. The number of onward infections
#' caused by each individual is drawn from a negative binomial distribution
#' parameterised by the mean R and dispersion parameter k.
#'
#' @param Rparams A list object containing parameters necessary for simulation
#' @param dt Size of the simulation time step
#' @param num_dt Number of simulation time steps. Simulation may end before this if the epidemic dies out.
#' @param min_epi_size Minimum number of infected individuals in a simulated outbreak.
#' @param maximum_attempts Number of attempts to simulate an outbreak greater than min_epi_size.
#' @param track_transmissions Whether or not to track who infected whom
#' @export
simulate_sir <- function(Rparams, dt, max_num_dt, min_epi_size, maximum_attempts, track_transmissions) {
    .Call('_EpiGenR_simulate_sir', PACKAGE = 'EpiGenR', Rparams, dt, max_num_dt, min_epi_size, maximum_attempts, track_transmissions)
}

#' Generate a transmission tree for an epidemic
#'
#' @param Repidemic A matrix with 4 columns: parent, time of infection, and time of recovery (e.g. output from simulate_sir)
#' @export
get_transmission_tree <- function(Repidemic) {
    .Call('_EpiGenR_get_transmission_tree', PACKAGE = 'EpiGenR', Repidemic)
}

