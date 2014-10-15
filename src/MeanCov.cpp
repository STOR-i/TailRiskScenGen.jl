#include <stdexcept>
#include <cassert>

#include "MeanCov.hpp"

const bool dbg = false;

MeanCov::MeanCov(arma::vec mean, arma::mat cov)
    : mean(mean), cov(cov), dim(mean.n_elem)
{
  if (mean.n_elem != cov.n_rows)
    throw std::logic_error("MeanCov: mean and covariance matrix "
			   " have incompatible dimensions.");
  if (!cov.is_square())
    throw std::logic_error("MeanCov: covariance must be a square "
			   "matrix.");
}

MeanCov::MeanCov(int dim, double *mean, double *cov)
  : mean(mean, dim, true), cov(cov, dim, dim, true)
{}

MeanCov::MeanCov(std::string normal_fn) {
  // Read Normal parameters
  std::ifstream normal_input(normal_fn, std::ios::in);
  std::string buffer;
  if (!normal_input)
    throw std::runtime_error("Cannot open file " + normal_fn);
  
  normal_input >> buffer;
  assert(buffer=="dim");
  normal_input >> dim;
  if (dbg) std::cerr << "dim = " << dim << std::endl;

  normal_input >> buffer;
  assert(buffer == "mean");
  mean = arma::vec(dim);
  for (unsigned int i = 0; i < dim; ++i) normal_input >> mean[i];
  if (dbg) std::cerr << mean << std::endl;

  normal_input >> buffer;
  assert(buffer == "cov");
  cov = arma::mat(dim, dim);
  for (unsigned int i = 0; i < dim; ++i) 
    for (unsigned int j = 0; j < dim; ++j) normal_input >> cov(i,j);
  normal_input.close();
}

void MeanCov::pprint (std::ostream &output) const {
  output << "dim " << dim;
  output << "\nmean " << mean.t();
  output << "cov\n" << cov;
}
