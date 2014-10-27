#include <iostream>
#include <armadillo>
#include "MVNormal.hpp"
#include "AggregationSampling.hpp"

int main() {
  arma::vec mean {0.1, -0.05, 0.01};
  arma::mat Sigma = {0.5, 0.1, 0.0,
		     0.1, 0.7, 0.3,
		     0.0, 0.3, 0.5};
  Sigma.reshape(3,3);
  int num_gen = 3;
  arma::mat cone_A = arma::eye(3,3);
  MeanCov norm(mean, Sigma);
  MVNormal norm_sampler(norm);
  double alpha = 1.6448;     //  = \Phi^{-1} (0.95)
  int num_scen = 30;

  std::cout << "Testing armadillo constructor..." << std::endl;
  AggregationSampling scen_gen_1(norm, cone_A, alpha);
  ScenarioSet scen_1 = scen_gen_1(num_scen);
  scen_1.pprint();

  std::cout << "Testing function constructor..." << std::endl;
  MVNormal normal_sampler(norm);
  AggregationSampling scen_gen_2(norm_sampler, norm, cone_A.memptr(), num_gen, alpha);
  ScenarioSet scen_2 = scen_gen_2(num_scen);
  scen_2.pprint();
  
  std::cout << "Testing copy constructor..." << std::endl;
  AggregationSampling scen_gen_3(scen_gen_1);
  ScenarioSet scen_3 = scen_gen_3(num_scen);
  scen_3.pprint();
  
  std::cout << "Testing scenario aggregation..." << std::endl;
  ScenarioSet scen_4 = rnorm_ScenarioSet(norm, num_scen);
  scen_4.pprint();
  ScenarioSet scen_5 = scen_gen_1(scen_4);
  scen_5.pprint();
}
