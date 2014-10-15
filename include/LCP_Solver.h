/* Siconos-Numerics, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#ifndef LCP_SOLVER_H
#define LCP_SOLVER_H

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

/** Structure used to handle with matrix in Numerics (interface to double*, SparseBlockStructuredMatrix and so on) \n
    Warning: one and only one storage is allowed and thus only one of the pointers below can be different from NULL
    \param storageType, int that identifies the type of storage (0: double*, 1:SparseBlockStructuredMatrix)
    \param size0, number of rows
    \param size1, number of columns
    \param matrix0 dense storage
    Related functions: prod(), subRowProd(), freeNumericsMatrix(), display()
*/

  typedef struct
  {
    int size0;
    int size1;
    double* matrix0;
  } NumericsMatrix;


  /** Matrix - vector product y = alpha*A*x + beta*y
      \param[in] sizeX dim of the vector x
      \param[in] sizeY dim of the vector y
      \param[in] alpha coefficient
      \param[in] A the matrix to be multiplied
      \param[in] x the vector to be multiplied
      \param[in] beta coefficient
      \param[in,out] y the resulting vector
  */
  void prodNumericsMatrix(int sizeX, int sizeY, double alpha, const NumericsMatrix* const A, 
			  const double* const x, double beta, double* y);



  typedef struct _SolverOptions
  {
    int isSet;                               /**< isSet int equal to false(0) if the parameters below have not been set (ie need to read default values) else true(1)*/
    int iSize;                               /**< iSize size of vector iparam */
    int * iparam;                            /**< iparam a list of int parameters (depends on each solver, see solver doc)*/
    int dSize;                               /**< dSize size of vector dparam */
    double * dparam;                         /**< dparam a list of double parameters (depends on each solver, see solver doc)*/
    int filterOn;                            /**< filterOn 1 to check solution validity after the driver call, else 0. Default = 1. (For example if
					      * filterOn = 1 for a LCP, lcp_compute_error() will be called at the end of the process) */
    int verboseMode;       /**< numericsOptions global options for numerics (verbose mode ...)*/
  } SolverOptions;
  
  /** delete the solver parameters :
      delete iparam and dparam;
      \param options the structure to be destroyed
  */
  void deleteSolverOptions(SolverOptions * options);


  /** \struct LinearComplementarityProblem LinearComplementarityProblem.h
   *  \brief Structure that contains and defines  \ref LCProblem
   *
   *   Find \f$(z,w)\f$ such that:\n
   *   \f{equation*}{
   *   \begin{cases}
   *   M \ z + q = w \\
   *   0 \le w \perp z \ge 0 \\
   *   \end{cases}
   *   \f}
   *
   * where \f$ w, z, q\f$ are vectors of size \f$n\f$ and \f$ M \f$ is a \f$n\times n\f$ matrix.
   * See \ref LCProblem for more details.
   */
  typedef struct
  {
    int size; /**<  size of the problem */
    NumericsMatrix* M ;/**< M matrix of the LCP (see the mathematical description)*/
    double * q;/**< vector of the LCP (see the mathematical description)*/
  } LinearComplementarityProblem;


  /** \fn   int linearComplementarity_driver(LinearComplementarityProblem* problem, double *z , double *w, SolverOptions* options)
   *
   General interface to solvers for Linear Complementarity Problems
    \param[in] problem the LinearComplementarityProblem structure which handles the problem (M,q)
    \param[in,out] z a n-vector of doubles which contains the solution of the problem.
    \param[in,out] w a n-vector of doubles which contains the solution of the problem.
    \param[in]  global_options the global options of Numerics
    \return info termination value
    - 0 : successful\n
    - >0 : otherwise see each solver for more information about the log info
    \author Franck Perignon
  */
  int linearComplementarity_driver(LinearComplementarityProblem* problem, double *z , double *w, SolverOptions* options);


  /** This function computes the input vector \f$ w = Mz + q \f$ and checks the validity of the vector z as a solution \n
   * of the LCP : \n
   * \f$
   *    0 \le z \perp Mz + q \ge 0
   * \f$
   * The criterion is based on \f$ \sum [ (z[i]*(Mz+q)[i])_{pos} + (z[i])_{neg} + (Mz+q)[i])_{neg} ] \f$ \n
   * with \f$ x_{pos} = max(0,x) \f$ and \f$ xneg = max(0,-x)\f$. \n
   * This sum is divided by \f$ \|q\| \f$ and then compared to tol.\n
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[in] tolerance threshold used to validate the solution: if the error is less than this value, the solution is accepted
   * \param[out] error the actual error of the solution with respect to the problem
   * \param[in] verbose print error if > than tolerance
   * \return status: 0 : convergence, 1: error > tolerance
   * \author Pascal Denoyelle, Franck Perignon
   */
  int lcp_compute_error(LinearComplementarityProblem* problem, double *z , double *w, 
			double tolerance, double* error, int verbose);


  /** lcp_lexicolemke is a direct solver for LCP based on pivoting method principle for degenerate problem \n
   * Choice of pivot variable is performed via lexicographic ordering \n
   *  Ref: "The Linear Complementarity Problem" Cottle, Pang, Stone (1992)\n
   * \param[in] problem structure that represents the LCP (M, q...)
   * \param[in,out] z a n-vector of doubles which contains the initial solution and returns the solution of the problem.
   * \param[in,out] w a n-vector of doubles which returns the solution of the problem.
   * \param[out] info an integer which returns the termination value:\n
   * 0 : convergence\n
   * 1 : iter = itermax\n
   * 2 : negative diagonal term\n
   * \param[in,out] options structure used to define the solver and its parameters.
   *
   *\author Mathieu Renouf
   */
  void lcp_lexicolemke(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LinearComplementarity
      \param options the pointer to the array of options to set
  */
  int linearComplementarity_lexicolemke_setDefaultSolverOptions(SolverOptions* options);

  /** screen display of solver parameters
      \param options the structure to be displayed
  */
  void printSolverOptions(SolverOptions* options);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
