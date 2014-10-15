#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include "ScenarioSet.hpp"
#include "MVNormal.hpp"

ScenarioSet::ScenarioSet(arma::mat scenarios, arma::vec probs)
  : scenarios(scenarios), probs(probs),
    dim(scenarios.n_rows), num_scen(scenarios.n_cols)
{
  if (scenarios.n_cols != probs.n_elem)
    throw std::logic_error("Scenario matrix must have number of"
			   " columns equal to the dimension of"
			   " probability vector.");
}

ScenarioSet::ScenarioSet(arma::mat scenarios)
  : scenarios(scenarios), probs(scenarios.n_cols),
    dim(scenarios.n_rows), num_scen(scenarios.n_cols)
{
  for (int s = 0; s < num_scen; ++s)
    probs[s] = 1.0/num_scen;
}

ScenarioSet::ScenarioSet(double *scenarios, double *probs, int dim, int num_scen)
  : scenarios(scenarios, dim, num_scen, true), probs(probs, num_scen, true), dim(dim), num_scen(num_scen)
{}

void ScenarioSet::pprint(std::ostream &output) const {
  const int width = 13;
  output.precision(3);
  output << std::left << std::setw(width) << "";
  for (int i = 0; i < dim; ++i)
    output << std::setw(width) << "Dim_" + std::to_string(i);
  output << std::setw(width) << "Prob." << std::endl;
  for (int s = 0; s < num_scen; ++s) {
    output << std::setw(width) << s;
    for (int i = 0; i < dim; ++i)
      output << std::setw(width) << std::scientific << get_scenario(s,i);
    output << std::setw(width) << std::scientific << 
      get_probability(s) << std::endl;
  }
}

void ScenarioSet::read(std::string fn) {
  std::ifstream input(fn, std::ios::in);
  std::string buffer;
  if (!input)
    throw std::runtime_error("Cannot open file " + fn);

  // Count lines to get num_scen
  num_scen = std::count(std::istreambuf_iterator<char>(input), 
			std::istreambuf_iterator<char>(), '\n') - 1;
  input.seekg(0);

  // Work out dimension of scenario set
  int counter = 0;
  while (input >> buffer) {
    if (buffer == "0") break;
    ++counter;
  }
  input.seekg(0);
  dim = counter - 1;

  std::getline(input, buffer);
  scenarios = arma::mat(dim, num_scen);
  probs = arma::vec(num_scen);
  for (int s = 0; s < num_scen; ++s) {
    input >> buffer;
    for (int i = 0; i < dim; ++i) input >> scenarios(i,s);
    input >> probs(s);
  }
}

ScenarioSet rnorm_ScenarioSet(const MeanCov &norm, int num_scen)
{
  MVNormal rng(norm);
  int dim = norm.get_dim();
  arma::mat scenarios(dim, num_scen);
  arma::vec probs(num_scen);

  double *ptr = scenarios.memptr();
  for (int i = 0; i < num_scen; ++i) {
    rng(ptr);
    ptr += dim;
    probs(i) = 1.0/num_scen;
  }
  return ScenarioSet(scenarios, probs);
}
