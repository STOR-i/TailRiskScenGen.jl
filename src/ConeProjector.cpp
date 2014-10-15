#include <cmath>
#include <iostream>
#include <iomanip>
#include <cblas.h>
#include <cstring>
#include "ConeProjector.hpp"

ConeProjector::ConeProjector(int num_gen, int dim, double * basis)
  : n(num_gen), m(dim),
    A(new double [m*n]), AtA(new double[n*n]),
    options(new SolverOptions)
{
  // Set Siconos options
  linearComplementarity_lexicolemke_setDefaultSolverOptions(options);
  // printSolverOptions(options);
  // Copy basis and calculate AtA
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      A[i*n + j] = basis[j*m + i];
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, 
	      n, m, 1.0, A, n, A, n, 0.0, AtA, n);
}

ConeProjector::ConeProjector(const ConeProjector &other)
  : n(other.n), m(other.m), A(new double[m*n]), AtA(new double[n*n]),
    options(new SolverOptions)
{
  // Set Siconos options
  linearComplementarity_lexicolemke_setDefaultSolverOptions(options);
  // Copy A and AtA
  memcpy(A, other.A, sizeof(double) * n * m);
  memcpy(AtA, other.AtA, sizeof(double) * n * n);
}
		
ConeProjector::ConeProjector(ConeProjector &&other)
  : n(other.n), m(other.m), A(other.A), AtA(other.AtA),
    options(other.options)
{
  other.A = NULL;
  other.AtA = NULL;
  other.n = 0;
  other.m = 0;
  other.options = NULL;
}

// Move assignment operator.
ConeProjector& ConeProjector::operator=(ConeProjector&& other)
{
   if (this != &other)
   {
      // Free the existing resource.
      delete[] A;
      delete[] AtA;
      delete options;
      deleteSolverOptions(options);

      // Copy the data pointer and its length from the 
      // source object.
      n = other.n;
      m = other.m;
      A = other.A;
      AtA = other.AtA;
      options = other.options;

      // Release the data pointer from the source object so that
      // the destructor does not free the memory multiple times.
      other.n = 0;
      other.m = 0;
      other.A = NULL;
      other.AtA = NULL;
      other.options = NULL;
   }
   return *this;
}

ConeProjector::ConeProjector()
  : n(0), m(0), A(NULL), AtA(NULL),
    options(NULL)
{}

ConeProjector::~ConeProjector()
{
  delete [] A;
  delete [] AtA;
  deleteSolverOptions(options);
  delete options;
}

void ConeProjector::project(const double *p, double *proj) {
  LinearComplementarityProblem prob;
  double z[n], w[n], q[n];
  cblas_dgemv(CblasRowMajor, CblasTrans,
	      m, n, -1.0, A, n, p, 1, 0.0, q, 1);
  NumericsMatrix mat;
  mat.size0 = n;
  mat.size1 = n;
  mat.matrix0 = AtA;
  prob.M = &mat;
  prob.q = q;
  prob.size = n;
  // linearComplementarity_display(&prob);
  linearComplementarity_driver(&prob, z, w, options);
  cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A, n,
	      z, 1, 0.0, proj, 1);
}

void ConeProjector::display()
{
  std::cout << "\n------------------------------------------" << std::endl;
  std::cout << "Number of generators: " << n << std::endl;
  std::cout << "Dimension: " << m << std::endl;
  int k = std::max(n,m);
  std::cout << "A" << std::setw(7*n - 1) << "\tAtA" << std::endl;

  for (int i = 0; i < k; ++i) {
    if (i < m) {
      for (int j = 0; j < n; ++j)
	std::cout << std::setprecision(2) << std::setw(6) << A[i*n + j];
    }
    else
      std::cout << std::setw(6 * n) << "";
    std::cout << "\t";

    if (i < n) {
      for (int j = 0; j < n; ++j)
	std::cout << std::setw(6) << AtA[i*n + j];
    }
    std::cout << std::endl;
  }
  std::cout << "------------------------------------------" << std::endl;
}
