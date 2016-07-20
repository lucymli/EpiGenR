#include <omp.h>
#include <Rcpp.h>
using namespace Rcpp;

//' Generate a transmission tree for an epidemic
//'
//' @param Repidemic A matrix with 4 columns: parent, time of infection, and time of recovery (e.g. output from simulate_sir)
//' @export
// [[Rcpp::export]]
NumericMatrix get_transmission_tree(SEXP Repidemic) {
  NumericMatrix epidemic(Repidemic);
  // Count how many edges there are in the transmission tree
  long nrow = epidemic.nrow();
  long counter = 0;
  for (int i=0; i!=epidemic.nrow(); ++i) {
    if (epidemic(i, 0) == -1) {
      ++counter;
      --nrow;
    }
  }
  NumericMatrix transmission_tree(nrow, 3);
  colnames(transmission_tree) = CharacterVector::create("from", "to", "length");
  double curr_infection_time;
  double parent_infection_time;
  double edge_length;
  int threads = omp_get_max_threads();
  omp_set_num_threads(threads);
#pragma omp parallel for
  for (int i=counter; i<epidemic.nrow(); ++i) {
    long parent_id = epidemic(i, 0);
    if (parent_id > -1) {
      curr_infection_time = epidemic(i, 1);
      parent_infection_time = epidemic(parent_id, 1);
      edge_length = curr_infection_time - parent_infection_time;
      transmission_tree(i-counter, 0) = parent_id;
      transmission_tree(i-counter, 1) = i;
      transmission_tree(i-counter, 2) = edge_length;
    }
  }
  return (transmission_tree);
}

/*
 ** R
 sourceCpp('../src/simulate.cpp')
 R0 <- 2
 k <- 0.5
 Tg <- 5
 N <- 5000
 S <- 4999
 dt <- 0.1
 total_dt <- 1000
 min_epi_size <- 20
 max_attempts <- 100
 params <- c(R0=R0, k=k, Tg=Tg, N=N, S=S, I=N-S)
 seed.num <- 101010133
 set.seed(seed.num)
 sim.outbreak <- simulate_sir(params, dt, total_dt, min_epi_size, max_attempts, TRUE)
 */
