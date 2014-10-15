#ifndef CONE_PROJECTOR_HPP
#define CONE_PROJECTOR_HPP

#include <LCP_Solver.h>

/**
 * Function object for projecting
 * points onto a finitely generated cone.
 */
class ConeProjector {
public:
  /**
   * Constructor
   * \param num_gen number of cone generators
   * \param dim dimension of space
   * \param basis dim \times num_gen matrix generating cone.
   *                         Stored in column major format.
   * \param opt solver options for underlying LCP
   */
  ConeProjector(int num_gen, int dim, double *basis);

  /**
   * Default constructor
   */
  ConeProjector();

  /**
   * Copy constructor
   */
  ConeProjector(const ConeProjector &);

  /**
   * Move constructor
   */
  ConeProjector(ConeProjector &&other);

  /**
   * Move assignment
   */
  ConeProjector &operator=(ConeProjector &&other);

  /**
   * Destructor
   */
  ~ConeProjector();
  
  /**
   * Projects point onto a finitely generated cone
   * \param[in] p point to be projected
   * \param[out] proj projection
   */
  void project(const double *p, double *proj);

  /**
   * Print projection data. Provided for
   * debugging purposes.
   */
  void display();

private:
  int n, m;
  double *A;
  double *AtA;
  SolverOptions *options;
};

#endif
