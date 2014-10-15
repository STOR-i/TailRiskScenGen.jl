#include <cstdio>
#include <cstring>
#include <cblas.h>
#include <iostream>
#include <stdexcept>
#include "ScenarioGenerationBase.hpp"
#include "ConeProjector.hpp"

ScenarioGenerationBase::ScenarioGenerationBase(const MeanCov &norm,
					       arma::mat cone_A, 
					       double alpha)
  : mean(norm.get_mean()), alpha(alpha)
{
  // Check arguments
  if (cone_A.n_rows != norm.get_dim())
    throw std::logic_error("Number of rows in cone basis must equal"
  			   " the dimension of the problem");
  dim = mean.n_elem;
  num_gen = cone_A.n_cols;
  P = arma::chol(norm.get_cov());
  inv_P = P.i();
  arma::mat A = P * cone_A;
  proj = new ConeProjector(num_gen, dim, A.memptr());
}

ScenarioGenerationBase::ScenarioGenerationBase(const MeanCov &norm,
					       double *cone_A, int num_gen,
					       double alpha)
  : num_gen(num_gen), mean(norm.get_mean()), alpha(alpha)
{
  dim = mean.n_elem;
  arma::mat arma_cone_A(cone_A, dim, num_gen, false, true);
  P = arma::chol(norm.get_cov());
  inv_P = P.i();
  arma::mat A = P * arma_cone_A;
  proj = new ConeProjector(num_gen, dim, A.memptr());
}

ScenarioGenerationBase::ScenarioGenerationBase(const ScenarioGenerationBase &other) 
  : dim(other.dim), num_gen(other.num_gen), mean(other.mean),
    P(other.P), inv_P(other.inv_P), 
    proj(new ConeProjector(*(other.proj))), alpha(other.alpha)
{}

ScenarioGenerationBase::~ScenarioGenerationBase() {
  delete proj;
}

void ScenarioGenerationBase::pprint(std::ostream &output) const {
  printf("Dimension: %d\n", this->dim);
  printf("Number of generators: %d\n", this->num_gen);
  output << mean << std::endl;
  output << P << std::endl;
  proj->display();
}

double ScenarioGenerationBase::projected_distance(arma::vec y) const {
  double *y_proj = new double[dim];
  arma::vec y_trans = inv_P.t() * (y - mean);
  proj->project(y_trans.memptr(), y_proj);
  double len = cblas_dnrm2(dim, y_proj, 1);
  delete [] y_proj;
  return len;
}

bool ScenarioGenerationBase::in_agg_region(arma::vec y) const {
  double len = projected_distance(y);
  if (len > alpha)
    return false;
  return true;
}
