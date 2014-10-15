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
  arma::mat cone_A = arma::eye(3,3);
  MeanCov norm(mean, Sigma);
  double beta = 0.95;
  int num_scen = 30;
  
  AggregationSampling scen_gen1(norm, cone_A, beta);
  AggregationSampling scen_gen(scen_gen1);
  ScenarioSet scen = scen_gen(num_scen);
  scen.pprint();
  ScenarioSet scen2 = scen_gen(scen);
  scen2.pprint();
  
  std::cout << "Agg sampling successful" << std::endl;
  MVNormal sampler(norm);
  ScenarioSet ss = rnorm_ScenarioSet(norm, 10);
  ss.pprint();
  ScenarioSet scen3 = scen_gen(ss);
  scen3.pprint();
}
