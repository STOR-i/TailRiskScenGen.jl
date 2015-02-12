#include "LCP_Solver.h"

/**
 * Find w >= 0, z >= 0 in R^n such that
 * M z + q = w
 * where M in R^(n x n) and q in R^n
 *
 * \param[out] w array to store solution for w
 * \param[out] z array to store solution for z
 * \param[in] M pointer to (n x n)-dimensional array
 * \param[in] q pointer to n-dimension array
 * \param[in] n dimension of problem
 * \return 1 for success and 0 for failure
 */
int lcp_julia_solve(double *w, double *z, const double *M, const double *q, int n) {
  int return_code;
  NumericsMatrix mat;
  LinearComplementarityProblem prob;
  SolverOptions options;
  linearComplementarity_lexicolemke_setDefaultSolverOptions(&options);
  mat.size0 = n;
  mat.size1 = n;
  mat.matrix0 = M;
  prob.M = &mat;
  prob.q = q;
  prob.size = n;
  return_code = linearComplementarity_driver(&prob, z, w, &options);
  deleteSolverOptions(&options);  
  return return_code;
}
