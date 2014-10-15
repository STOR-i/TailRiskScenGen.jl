#include <iostream>
#include <fstream>
#include <armadillo>
#include "MVNormal.hpp"

int main() {
  const int num_scen = 10000;
  const int dim = 3;
  arma::vec mean = {0.0, 0.0, 0.0};
  arma::mat cov = {1.0, 0.1, 0.0,
		   0.1, 1.0, 0.0,
		   0.0, 0.0, 2.0};
  cov.reshape(3,3);
  MeanCov norm_param(mean, cov);
  std::ofstream out("norm.txt");
  norm_param.pprint();

  MVNormal norm(norm_param);
  
  arma::mat::fixed<dim, num_scen> scenarios;
  double *ptr = scenarios.memptr();

  for (int s = 0; s < num_scen; ++s) {
    norm(ptr);
    ptr += dim;
  }

  std::cout << arma::cov(scenarios.t());
}
