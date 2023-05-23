#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List compute_gradients_bptt(const arma::mat& a_x, const arma::mat& y_feedback_activation, const arma::mat& x_feedback_activation, const List& weights, double c) {
  arma::mat w_xh = as<arma::mat>(weights["w_xh"]);
  arma::mat w_hh = as<arma::mat>(weights["w_hh"]);
  arma::mat w_ho = as<arma::mat>(weights["w_ho"]);

  // Initialize the gradients
  arma::mat dw_xh = arma::zeros<arma::mat>(w_xh.n_rows, w_xh.n_cols);
  arma::mat dw_hh = arma::zeros<arma::mat>(w_hh.n_rows, w_hh.n_cols);
  arma::mat dw_ho = arma::zeros<arma::mat>(w_ho.n_rows, w_ho.n_cols);

  // Implement the BPTT algorithm to compute the gradients
  // Assuming a single time step
  arma::mat delta_out = y_feedback_activation - x_feedback_activation;
  arma::mat delta_hidden = w_ho.t() * delta_out % a_x % (1 - a_x); // Element-wise multiplication (%)

  dw_xh += delta_hidden * a_x.t();
  dw_hh += delta_hidden * a_x.t(); // Assuming no previous hidden state
  dw_ho += delta_out * a_x.t();

  return List::create(
    _["dw_xh"] = dw_xh,
    _["dw_hh"] = dw_hh,
    _["dw_ho"] = dw_ho
  );
}

