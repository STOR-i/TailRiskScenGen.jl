#include <stdlib.h>
#include <stddef.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <cblas.h>
#include "LCP_Solver.h"
#include "debug.h"

void numericsError(char * functionName, char* message)
{
  char output[200] = "Numerics error - ";
  strcat(output, functionName);
  strcat(output, message);
  strcat(output, ".\n");
  fprintf(stderr, "%s", output);
  exit(EXIT_FAILURE);
}

void prodNumericsMatrix(int sizeX, int sizeY, double alpha, const NumericsMatrix* const A, const double* const x, double beta, double* y)
{
  assert(A);
  assert(x);
  assert(y);
  assert(A->size0 == sizeY);
  assert(A->size1 == sizeX);
  cblas_dgemv(CblasColMajor, CblasNoTrans, sizeY, sizeX, alpha, A->matrix0, sizeY, x, 1, beta, y, 1);
}

void deleteSolverOptions(SolverOptions* op)
{
  if(op)
  {
    if (op->iparam != NULL)
      free(op->iparam);
    op->iparam = NULL;
    if (op->dparam != NULL)
      free(op->dparam);
    op->dparam = NULL;
  }
}

void lcp_lexicolemke(LinearComplementarityProblem* problem, double *zlem , double *wlem , int *info , SolverOptions* options)
{
  /* matrix M of the lcp */
  double * M = problem->M->matrix0;
  assert(M);
  /* size of the LCP */
  int dim = problem->size;
  assert(dim>0);
  int dim2 = 2 * (dim + 1);

  int i, drive, block, Ifound;
  int ic, jc;
  int ITER;
  int nobasis;
  int itermax = options->iparam[0];

  i=0;
  int n = problem->size;
  double *q = problem->q;
  
  while ((i < (n - 1)) && (q[i] >= 0.)) 
    i++;
  
  if ((i == (n - 1)) && (q[n - 1] >= 0.))
  {
    /* TRIVIAL CASE : q >= 0
     * z = 0 and w = q is solution of LCP(q,M)
     */
    for (int j = 0 ; j < n; j++)
    {
      zlem[j] = 0.0;
      wlem[j] = q[j];
    }
    *info = 0;
    options->iparam[1] = 0;   /* Number of iterations done */
    options->dparam[1] = 0.0; /* Error */
    if (options->verboseMode > 0)
      printf("lcp_lexicolemke: found trivial solution for the LCP (positive vector q => z = 0 and w = q). \n");
    return ;
  }
  
  double z0, zb, dblock;
  double pivot, tovip;
  double tmp;
  int *basis;
  double** A;

  /*output*/
  options->iparam[1] = 0;

  /* Allocation */
  basis = (int *)malloc(dim * sizeof(int));
  A = (double **)malloc(dim * sizeof(double*));

  for (ic = 0 ; ic < dim; ++ic)
    A[ic] = (double *)malloc(dim2 * sizeof(double));

  /* construction of A matrix such that
   * A = [ q | Id | -d | -M ] with d = (1,...1)
   */
  /* We need to init only the part corresponding to Id */
  for (ic = 0 ; ic < dim; ++ic)
    for (jc = 1 ; jc <= dim; ++jc)
      A[ic][jc] = 0.0;

  for (ic = 0 ; ic < dim; ++ic)
    for (jc = 0 ; jc < dim; ++jc)
      A[ic][jc + dim + 2] = -M[dim * jc + ic];

  assert(problem->q);

  for (ic = 0 ; ic < dim; ++ic) A[ic][0] = problem->q[ic];

  for (ic = 0 ; ic < dim; ++ic) A[ic][ic + 1 ] =  1.0;
  for (ic = 0 ; ic < dim; ++ic) A[ic][dim + 1] = -1.0;

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("%1.2e ", A[i][j]) }
      DEBUG_PRINT("\n")});
  /* End of construction of A */

  Ifound = 0;


  for (ic = 0 ; ic < dim  ; ++ic) basis[ic] = ic + 1;

  drive = dim + 1;
  block = 0;
  z0 = A[block][0];
  ITER = 0;

  /* Start research of argmin lexico */
  /* With this first step the covering vector enter in the basis */
  for (ic = 1 ; ic < dim ; ++ic)
  {
    zb = A[ic][0];
    if (zb < z0)
    {
      z0    = zb;
      block = ic;
    }
    else if (zb == z0)
    {
      for (jc = 0 ; jc < dim ; ++jc)
      {
        dblock = A[block][1 + jc] - A[ic][1 + jc];
        if (dblock < 0)
        {
          break;
        }
        else if (dblock > 0)
        {
          block = ic;
          break;
        }
      }
    }
  }

  /* Stop research of argmin lexico */
  DEBUG_PRINTF("Pivoting %i and %i\n", block, drive);

  pivot = A[block][drive];
  tovip = 1.0 / pivot;

  /* Pivot < block , drive > */
  A[block][drive] = 1;
  for (ic = 0       ; ic < drive ; ++ic) A[block][ic] = A[block][ic] * tovip;
  for (ic = drive + 1 ; ic < dim2  ; ++ic) A[block][ic] = A[block][ic] * tovip;

  /* */

  for (ic = 0 ; ic < block ; ++ic)
  {
    tmp = A[ic][drive];
    for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
  }
  for (ic = block + 1 ; ic < dim ; ++ic)
  {
    tmp = A[ic][drive];
    for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
  }

   nobasis = basis[block];
  basis[block] = drive;

  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));
  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("%1.2e ", A[i][j]) }
      DEBUG_PRINT("\n")});

  while (ITER < itermax && !Ifound)
  {

    ++ITER;

    if (nobasis < dim + 1)      drive = nobasis + (dim + 1);
    else if (nobasis > dim + 1) drive = nobasis - (dim + 1);

    DEBUG_PRINTF("driving variable %i \n", drive);

    /* Start research of argmin lexico for minimum ratio test */
    pivot = 1e20;
    block = -1;

    for (ic = 0 ; ic < dim ; ++ic)
    {
      zb = A[ic][drive];
      if (zb > 0.0)
      {
        z0 = A[ic][0] / zb;
        if (z0 > pivot) continue;
        if (z0 < pivot)
        {
          pivot = z0;
          block = ic;
        }
        else
        {
          for (jc = 1 ; jc < dim + 1 ; ++jc)
          {
            assert(block >=0 && "lcp_lexicolemke: block <0");
            dblock = A[block][jc] / pivot - A[ic][jc] / zb;
            if (dblock < 0.0) break;
            else if (dblock > 0.0)
            {
              block = ic;
              break;
            }
          }
        }
      }
    }
    if (block == -1)
    {
      Ifound = 1;
      DEBUG_PRINT("The pivot column is nonpositive !\n"
          "It either means that the algorithm failed or that the LCP is infeasible\n"
          "Check the class of the M matrix to find out the meaning of this\n");
      break;
    }

    if (basis[block] == dim + 1) Ifound = 1;

    /* Pivot < block , drive > */
    pivot = A[block][drive];
    tovip = 1.0 / pivot;
    A[block][drive] = 1;

    for (ic = 0       ; ic < drive ; ++ic) A[block][ic] = A[block][ic] * tovip;
    for (ic = drive + 1 ; ic < dim2  ; ++ic) A[block][ic] = A[block][ic] * tovip;

    /* */

    for (ic = 0 ; ic < block ; ++ic)
    {
      tmp = A[ic][drive];
      for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
    }
    for (ic = block + 1 ; ic < dim ; ++ic)
    {
      tmp = A[ic][drive];
      for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
    }

    nobasis = basis[block];
    basis[block] = drive;

    DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

    DEBUG_PRINT("total matrix\n");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("%1.2e ", A[i][j]) }
      DEBUG_PRINT("\n")});

  } /* end while*/

  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("%1.2e ", A[i][j]) }
      DEBUG_PRINT("\n")});

  for (ic = 0 ; ic < dim; ++ic)
  {
    drive = basis[ic];
    if (drive < dim + 1)
    {
      zlem[drive - 1] = 0.0;
      wlem[drive - 1] = A[ic][0];
    }
    else if (drive > dim + 1)
    {
      zlem[drive - dim - 2] = A[ic][0];
      wlem[drive - dim - 2] = 0.0;
    }
  }

  options->iparam[1] = ITER;

  if (Ifound) *info = 0;
  else *info = 1;

  free(basis);

  for (i = 0 ; i < dim ; ++i) free(A[i]);
  free(A);
}

int linearComplementarity_driver(LinearComplementarityProblem* problem, double *z , double *w, SolverOptions* options)
{
  /********************
   * 0 - Check inputs *
   ********************/

  if (options == NULL)
    numericsError("lcp_driver", "null input for solver options");
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("lcp_driver", "null input for LinearComplementarityProblem and/or unknowns (z,w)");

  int NoDefaultOptions = options->isSet; /* true(1) if the SolverOptions structure has been filled in else false(0) */

  if (NoDefaultOptions == 0)
  {
    numericsError("lcp_driver_DenseMatrix", "options for solver have not been set");
  }

  if (options->verboseMode > 0)
    printSolverOptions(options);

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;

  /******************************************
   *  1 - Check for trivial solution
   ******************************************/

  int i = 0;
  int n = problem->size;
  double *q = problem->q;
/*  if (!((options->solverId == SICONOS_LCP_ENUM) && (options->iparam[0] == 1 )))*/
    {      
      while ((i < (n - 1)) && (q[i] >= 0.)) i++;
      if ((i == (n - 1)) && (q[n - 1] >= 0.))
      {
        /* TRIVIAL CASE : q >= 0
         * z = 0 and w = q is solution of LCP(q,M)
         */
        for (int j = 0 ; j < n; j++)
        {
          z[j] = 0.0;
          w[j] = q[j];
        }
        info = 0;
        options->dparam[1] = 0.0; /* Error */
        if (options->verboseMode > 0)
          printf("LCP_driver_DenseMatrix: found trivial solution for the LCP (positive vector q => z = 0 and w = q). \n");
        return info;
      }
    }

  /*************************************************
   *  2 - Call Lemke solver (if no trivial sol.)
   *************************************************/

  if (options->verboseMode == 1)
    printf(" ========================== Call Lemke solver for Linear Complementarity problem ==========================\n");

  /****** Lemke algorithm ******/
  /* IN: itermax
     OUT: iter */
   lcp_lexicolemke(problem, z , w , &info , options);

  /*************************************************
   *  3 - Computes w = Mz + q and checks validity
   *************************************************/
  if (options->filterOn > 0)
  {
    int info_ = lcp_compute_error(problem, z, w, options->dparam[0], &(options->dparam[1]),
				  options->verboseMode);
    if (info <= 0) /* info was not setor the solver was happy */
      info = info_;
  }
  
  return info;
}

void lcp_compute_error_only(unsigned int n, double *z , double *w, double * error)
{
  /* Checks complementarity */

  *error = 0.;
  double zi, wi;
  for (unsigned int i = 0 ; i < n ; i++)
  {
    zi = z[i];
    wi = w[i];
    if (zi < 0.0)
    {
      *error += -zi;
      if (wi < 0.0) *error += zi * wi;
    }
    if (wi < 0.0) *error += -wi;
    if ((zi > 0.0) && (wi > 0.0)) *error += zi * wi;
  }
}

int lcp_compute_error(LinearComplementarityProblem* problem, double *z , double *w, double tolerance, 
		      double * error, int verbose)
{
  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("lcp_compute_error", "null input for problem and/or z and/or w");

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  unsigned int n = problem->size;
  cblas_dcopy(n , problem->q , incx , w , incy);  // w <-q
  prodNumericsMatrix(n, n, 1.0, problem->M, z, 1.0, w);
  double normq = cblas_dnrm2(n , problem->q , incx);
  lcp_compute_error_only(n, z, w, error);
  *error = *error / (normq + 1.0); /* Need some comments on why this is needed */
  if (*error > tolerance)
  {
    if (verbose > 0) 
	printf(" Numerics - lcp_compute_error : error = %g > tolerance = %g.\n", *error, tolerance);
    return 1;
  }
  else
    return 0;
}

int linearComplementarity_lexicolemke_setDefaultSolverOptions(SolverOptions* options)
{
  if (options->verboseMode > 0)
  {
    printf("Set the Default SolverOptions for the Lemke Solver\n");
  }

  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->verboseMode = 0;
  options->dparam[0] = 1e-6;
  options->iparam[0] = 10000;
  return 0;
}

void printSolverOptions(SolverOptions* options)
{
  printf("\n ========== Numerics Non Smooth Solver parameters: \n");
  if (options->isSet == 0)
    printf("The solver parameters have not been set. \t options->isSet = %i \n", options->isSet);
  else
  {
    printf("The solver parameters below have  been set \t options->isSet = %i\n", options->isSet);
    printf("Name of the solver\t\t\t\t Lemke \n");
    if (options->iparam != NULL)
    {
      printf("int parameters \t\t\t\t\t options->iparam\n");
      printf("size of the int parameters\t\t\t options->iSize = %i\n", options->iSize);
      for (int i = 0; i < options->iSize; ++i)
        printf("\t\t\t\t\t\t options->iparam[%i] = %d\n", i, options->iparam[i]);
    }
    if (options->dparam != NULL)
    {
      printf("double parameters \t\t\t\t options->dparam\n");
      printf("size of the double parameters\t\t\t options->dSize = %i\n", options->dSize);
      for (int i = 0; i < options->iSize; ++i)
        printf("\t\t\t\t\t\t options->dparam[%i] = %.6le\n", i, options->dparam[i]);
    }
  }

  printf("See Lemke documentation for parameters definition)\n");
  printf("\n");
}
