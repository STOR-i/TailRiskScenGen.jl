#ifndef MEANCOV_HPP
#define MEANCOV_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <armadillo>

/**
 * Struct representing a set of multivariate
 * Normal distribution parameters.
 */

struct MeanCov {
public:
  /**
   * Constructs a MeanCov object with
   * specified mean vector and covariance matrix.
   *
   * Construction includes checks for consistency
   * between mean and covariance.
   */
  MeanCov(arma::vec mean, arma::mat cov);

  /**
   * Constructs a MeanCov object with
   * specified mean vector and covariance matrix.
   *
   * Construction includes checks for consistency
   * between mean and covariance.
   */
  MeanCov(int dim, double *mean, double *cov);

  /**
   * Constructs MeanCov object by reading
   * a specified file.
   *
   * \param normal_fn name of file
   */
  MeanCov(std::string normal_fn);

  arma::vec get_mean() const {return mean;}

  arma::mat get_cov() const {return cov;}

  size_t get_dim() const {return dim;}

  /**
   * Prints string representation of MeanCov object
   * to a specified output stream
   */
  void pprint (std::ostream &out = std::cout) const;
private:
  arma::vec mean;
  arma::mat cov;
  size_t dim;
};

#endif
