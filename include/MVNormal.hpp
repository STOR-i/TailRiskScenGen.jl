#ifndef MVNORMAL_HPP
#define MVNORMAL_HPP

#include <random>
#include <chrono>
#include <armadillo>
#include "MeanCov.hpp"

using namespace std::chrono;

/**
 * Function object for sampling from a
 * multivariate Normal distribution
 */
class MVNormal {
public:

  /**
   * Constructor
   *
   * Random number generator is currently set based on time.
   * \param norm parameters of normal distribution
   */
  MVNormal(const MeanCov &norm)
    : mean(norm.get_mean())
  {
    P = arma::chol(norm.get_cov()).t();
    dim = norm.get_dim();
    system_clock::duration dur = high_resolution_clock::now().time_since_epoch();
    engine.seed(dur.count());
  }

  /**
   * Copy constructor
   */
  MVNormal(const MVNormal &other)
    : mean(other.mean), 
      P(other.P), dim(other.dim)
  {
    system_clock::duration dur = high_resolution_clock::now().time_since_epoch();
    engine.seed(dur.count());
  }

  /**
   * Sample a point from the specified multivariate Normal
   * distribution.
   */
  void operator()(double *vec)
  {
    arma::vec point(vec, dim, false, true);
    for (int i = 0; i < dim; ++i)
      point[i] = distr(engine);
    point = P * point + mean;
  }

private:
  std::mt19937 engine;
  std::normal_distribution<double> distr;
  arma::vec mean;
  arma::mat P;
  int dim;
};

#endif
