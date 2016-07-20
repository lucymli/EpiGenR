#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

//' Simulate an epidemic according to an SIR model. The number of onward infections
//' caused by each individual is drawn from a negative binomial distribution
//' parameterised by the mean R and dispersion parameter k.
//'
//' @param Rparams A list object containing parameters necessary for simulation
//' @param dt Size of the simulation time step
//' @param num_dt Number of simulation time steps. Simulation may end before this if the epidemic dies out.
//' @param min_epi_size Minimum number of infected individuals in a simulated outbreak.
//' @param maximum_attempts Number of attempts to simulate an outbreak greater than min_epi_size.
//' @param track_transmissions Whether or not to track who infected whom
//' @export
// [[Rcpp::export]]
List simulate_sir(SEXP Rparams, double dt, int max_num_dt, int min_epi_size,
                  int maximum_attempts, bool track_transmissions) {
  NumericVector params(Rparams);
  double R0 = params["R0"]; // Reproductive number at start of epidemic
  double Tg = params["Tg"]; // Mean duration of infectiousness (generation time)
  double k = params["k"]; // Dispersion parameter of a negative binomial offspring distribution
  int N = params["N"]; // Total population size
  int S = params["S"]; // Number of susceptible individuals at start of epidemic
  int I = params["I"]; // Number of infectious individuals at start of epidemic
  int total_infected = I; // Cumulative number of infected individuals
  NumericVector prevalence(max_num_dt, 0.0);
  prevalence[0] = I;
  NumericVector parents(N, 0.0);
  NumericVector infection_times(N, 0.0);
  NumericVector recovery_times(N, 0.0);
  LogicalVector currently_infected(N, false);
  for (int i=0; i<I; ++i) {
    parents[i] = -1;
    currently_infected[i] = true;
  }
  int num_dt_taken = 0;
  int new_infections = 0;
  int individual_R = 0;
  int num_recover = 0;
  for (; num_dt_taken < max_num_dt; ++num_dt_taken) {
    if (track_transmissions) {
      for (int j=0; j<total_infected; ++j) { //j = infector
        if (currently_infected[j]) {
          num_recover = rbinom(1, 1, dt/Tg)[0];
          if (num_recover > 0) {
            prevalence[num_dt_taken] -= num_recover;
            recovery_times[j] = num_dt_taken;
            individual_R = rnbinom_mu(k, R0*(double)S/(double)N);
            if ((individual_R < 1) & (total_infected < min_epi_size)) {
              int counter = maximum_attempts;
              while ((counter > 0) & (individual_R < 0)) {
                individual_R = rnbinom_mu(k, R0*(double)S/(double)N);
                --counter;
              }
            }
            individual_R = std::min(individual_R, S);
            if (individual_R > 0) {
              S -= individual_R;
              prevalence[num_dt_taken] += individual_R;
            }
            for (int m=total_infected; m < total_infected+individual_R; ++m) { // m = infected
              parents[m] = j;
              infection_times[m] = num_dt_taken;
              currently_infected[m] = true;
            }
            total_infected += individual_R;
            currently_infected[j] = false;
          }
        }
      }
    }
    else {
      num_recover = rbinom(1, prevalence[num_dt_taken], dt/Tg)[0];
      new_infections = rnbinom_mu(k*(double)num_recover, R0*(double)S/(double)N*(double)num_recover);
      if ((new_infections < 1) & (total_infected < min_epi_size)) {
        int counter = maximum_attempts;
        while ((counter > 0) & (new_infections < 0)) {
          new_infections = rnbinom_mu(k*(double)num_recover, R0*(double)S/(double)N*(double)num_recover);
          --counter;
        }
      }
      if (new_infections > S) new_infections = S;
      total_infected += new_infections;
      prevalence[num_dt_taken] += new_infections - num_recover;
      S -= new_infections;
    }
    if (num_dt_taken < (max_num_dt-1)) prevalence[num_dt_taken+1] = prevalence[num_dt_taken];
    if (prevalence[num_dt_taken] < 1) break;
    if (S < 1) break;
  }
  NumericVector prevalence_output (num_dt_taken, 0.0);
  for (int i=0; i!=num_dt_taken; ++i) prevalence_output[i] = prevalence[i];
  if (total_infected < min_epi_size) return(List::create(Named("prevalence")=-1));
  if (track_transmissions) {
    NumericMatrix output(total_infected, 3);
    colnames(output) = CharacterVector::create("parents", "infection_time", "recovery_time");
    for (int i=0; i!=total_infected; ++i) {
      output(i, 0) = parents[i];
      output(i, 1) = infection_times[i]*dt;
      if (currently_infected[i]) recovery_times[i] = num_dt_taken;
      output(i, 2) = recovery_times[i]*dt;
    }
    return List::create(Named("infected")=output,
                        Named("prevalence")=prevalence_output,
                        Named("params")=params,
                        Named("dt")=dt,
                        Named("total_dt")=num_dt_taken,
                        Named("total_infected")=total_infected);
  }
  return List::create(Named("prevalence")=prevalence_output,
                      Named("params")=params,
                      Named("dt")=dt,
                      Named("total_dt")=num_dt_taken,
                      Named("total_infected")=total_infected);
}


//' tmp <- simulate_sir(c(R0=5, Tg=5, k=1, N=1000, S=999, I=1), 0.1, 1000, 20, 1, FALSE)
