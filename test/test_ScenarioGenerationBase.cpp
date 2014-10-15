#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <string>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include "ScenarioGenerationBase.hpp"
#include "MeanCov.hpp"

const bool dbg=true;

void run_test(std::string fn);

int main() {
  run_test("test/results_1.txt");
  run_test("test/results_2.txt");
  run_test("test/results_3.txt");
  run_test("test/results_4.txt");
  run_test("test/results_5.txt");
  run_test("test/results_6.txt");
  return 0;
}


void run_test(std::string fn) {
  static int counter = 1;
  const int max_err_display = 5;
  int err_count = 0;

  double beta;
  int dim, num_gen, num_points;
  std::string buffer;

  std::ifstream input(fn, std::ios::in);
  if (!input) {
    std::string message = "Cannot open file " + fn + ". Try "
      "running test from project root.";
    throw std::runtime_error("Cannot open file " + fn);
  }

  input >> buffer;
  assert(buffer=="beta");
  input >> beta;
  // if (dbg) std::cout << "Beta = " << beta << std::endl;

  std::cout << "-----------\nTest " << counter << "\n"
	    << "-----------" << std::endl;
  std::cout << "Results file: " << fn << std::endl;

  input >> buffer;
  assert(buffer=="dim");
  input >> dim;
  if (dbg) std::cout << "dim = " << dim << std::endl;

  input >> buffer;
  assert(buffer=="num_gen");
  input >> num_gen;
  if (dbg) std::cout << "num_gen = " << num_gen << std::endl;

  input >> buffer;
  assert(buffer=="num_points");
  input >> num_points;
  if (dbg) std::cout << "num_points = " << num_points << std::endl;

  input >> buffer;
  assert(buffer == "mean");
  arma::vec mean(dim);
  for (int i = 0; i < dim; ++i) input >> mean[i];
  if (dbg) std::cout << mean << std::endl;

  input >> buffer;
  assert(buffer == "cov");
  arma::mat cov(dim, dim);
  for (int i = 0; i < dim; ++i) 
    for (int j = 0; j < dim; ++j) input >> cov(i,j);
  if (dbg) std::cout << cov << std::endl;

  input >> buffer;
  assert(buffer == "cone_basis");
  arma::mat cone_A(dim, num_gen);
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < num_gen; ++j)
      input >> cone_A[j*dim + i];
  if (dbg) std::cout << cone_A << std::endl;

  //MeanCov norm(mean, cov);
  MeanCov norm(dim, mean.memptr(), cov.memptr());
  //ScenarioGenerationBase sg(norm, cone_A, beta);
  ScenarioGenerationBase sg(norm, cone_A.memptr(), num_gen, beta);


  sg.pprint();
  double dist, true_dist;
  input >> buffer;
  assert(buffer == "points");
  arma::vec point(dim);
  for (int i = 0; i < num_points; ++i) {
    for (int j = 0; j < dim; ++j) input >> point[j];
    //if (dbg) std::cout << "Point: " << point << std::endl;
    input >> true_dist;
    dist = sg.projected_distance(point);
    if (fabs(true_dist - dist) > 1e-2) {
      if (err_count <= max_err_display) {
	std::cerr << "Failure at " << i << ": ";
	std::cerr << point << std::endl;
	std::cerr << "\nProjected distance: " << dist;
	std::cerr << "\nTrue projected distance: " << true_dist
		  << std::endl;
	std::cerr << "------------------------------------" << std::endl;
      }
      else if (err_count == max_err_display + 1)
	std::cerr << "Etc..." << std::endl;
      ++err_count;
    }
  }
  std::cout << num_points - err_count << "/" << num_points
	    << " points successful!\n" << std::endl;
  ++counter;
  
  input.close();
}
