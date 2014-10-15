#include <iostream>
#include <iomanip>
#include <cmath>
#include "ConeProjector.hpp"

bool vec_equals_vec(double *x, double *y, size_t dim, double eps);

void write_results_header(int dim);
void write_result(double *point, double *true_proj, double *res, int dim, double eps);
void test_1();
void test_2();
void test_3();

int main() {
  test_1();
  test_2();
  test_3();
  return 0;
}

void test_1() {
  std::cout << "\nTest 1..." << std::flush;
  double basis[4] = {0.0, 1.0, 1.0, 1.0};  // 0.0  1.0
                                           // 1.0  1.0
  double res[2];
  ConeProjector proj1(2, 2, basis);
  ConeProjector proj(proj1);
  proj.display();

  write_results_header(2);
  double y_1[2] = {-1, 3};
  double y_1_proj[2] = {0, 3};
  proj.project(y_1, res);
  write_result(y_1, y_1_proj, res, 2, 1e-5);
  
  double y_2[2] = {-1, -1};
  double y_2_proj[2] = {0, 0};
  proj.project(y_2, res);
  write_result(y_2, y_2_proj, res, 2, 1e-5);

  double y_3[2] = {2, 0};
  double y_3_proj[2] = {1, 1};
  proj.project(y_3, res);
  write_result(y_3, y_3_proj, res, 2, 1e-5);
  
  double y_4[2] = {0.25, 1.0};
  double y_4_proj[2] = {0.25, 1.0};
  proj.project(y_4, res);
  write_result(y_4, y_4_proj, res, 2, 1e-5);

}

void test_2() {
  std::cout << "\nTest 2..." << std::flush;
  double basis[6] = {0.0, 1.0, 1.0, 1.0, 0.5, 1.0};  // 0.0  1.0  0.5
                                                     // 1.0  1.0  1.0
  double res[2];
  ConeProjector proj1(3, 2, basis);
  ConeProjector proj(proj1);
  proj.display();

  write_results_header(2);
  double y_1[2] = {-1, 3};
  double y_1_proj[2] = {0, 3};
  proj.project(y_1, res);
  write_result(y_1, y_1_proj, res, 2, 1e-5);
  
  double y_2[2] = {-1, -1};
  double y_2_proj[2] = {0, 0};
  proj.project(y_2, res);
  write_result(y_2, y_2_proj, res, 2, 1e-5);

  double y_3[2] = {2, 0};
  double y_3_proj[2] = {1, 1};
  proj.project(y_3, res);
  write_result(y_3, y_3_proj, res, 2, 1e-5);
  
  double y_4[2] = {0.25, 1.0};
  double y_4_proj[2] = {0.25, 1.0};
  proj.project(y_4, res);
  write_result(y_4, y_4_proj, res, 2, 1e-5);
}

void test_3() {
  std::cout << "\nTest 3..." << std::flush;
  double basis[6] = {0.0, 0.0, 1.0, 1.0, 1.0, 1.0};  // 0.0  1.0
                                                     // 0.0  1.0
                                                     // 1.0  1.0
  double res[3];
  ConeProjector (2, 3, basis);
  ConeProjector proj;
  proj = ConeProjector(2, 3, basis);
  proj.display();

  write_results_header(3);
  double y_1[3] = {-1, -1, 1};
  double y_1_proj[3] = {0,0,1};
  proj.project(y_1, res);
  write_result(y_1, y_1_proj, res, 3, 1e-5);

  double y_2[3] = {3.0, 0.0, 0.0};
  double y_2_proj[3] = {1.0, 1.0, 1.0};
  proj.project(y_2, res);
  write_result(y_2, y_2_proj, res, 3, 1e-5);

  double y_3[3] = {-1.0, -1.0, -1.0};
  double y_3_proj[3] = {0.0, 0.0, 0.0};
  proj.project(y_3, res);
  write_result(y_3, y_3_proj, res, 3, 1e-5);

  double y_4[3] = {1.0, 2.0, 3.0};
  double y_4_proj[3] = {1.5, 1.5, 3.0};
  proj.project(y_4, res);
  write_result(y_4, y_4_proj, res, 3, 1e-5);

  double y_5[3] = {2.5, 2.5, 6.0};
  double y_5_proj[3] = {2.5, 2.5, 6.0};
  proj.project(y_5, res);
  write_result(y_5, y_5_proj, res, 3, 1e-5);
}


void write_results_header(int dim) {
  int col_width = 6 * dim;
  std::cout << "|" << std::setw(col_width) << "y"
	    << "|" << std::setw(col_width) << "True p_k(y)"
	    << "|" << std::setw(col_width) << "Calculation"
	    << "|" << std::endl;
  for (int i = 0; i < 4 * col_width + 4; ++i) std::cout << "-";
  std::cout << std::endl;
}


void write_result(double *point, double *true_proj, 
		  double *res, int dim, double eps) {
  std::cout << "|";
  for (int i = 0; i < dim; ++i)
    std::cout << std::setw(6) << std::setprecision(2) << std::fixed << point[i];
  std::cout << "|" ;
  for (int i = 0; i < dim; ++i)
    std::cout << std::setw(6) << true_proj[i];
  std::cout << "|" ;
  for (int i = 0; i < dim; ++i)
    std::cout << std::setw(6) << res[i];
  std::cout << "|" ;
  if (vec_equals_vec(res, true_proj, dim, eps))
    std::cout << std::setw(6 * dim) << "...success" << std::endl;
  else std::cout << std::setw(6 * dim) << "...fail" << std::endl;
}

bool vec_equals_vec(double *x, double *y, size_t dim, double eps) {
  double err = 0.0;
  for (unsigned int i = 0; i < dim; ++i)
    err += fabs(x[i] - y[i]);
  if (err > eps) return false;
  else return true;
}
