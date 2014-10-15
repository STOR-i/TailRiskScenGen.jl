#ifndef SCENARIO_GENERATION_BASE_HPP
#define SCENARIO_GENERATION_BASE_HPP

#include <iostream>
#include <armadillo>
#include "MeanCov.hpp"

class ConeProjector;

/**
 * Base class from which other scenario
 * aggregation methods are derived.
 *
 * Class currently assumes a Normal aggregation
 * region
 */
class ScenarioGenerationBase {
public:
  /**
   * Constructor
   *
   * \param norm object representing Normal parameters of returns
   *             distribution
   * \param cone_A matrix of cone generators for feasible region.
   *               Each column represents a generator
   * \param alpha quantile value (rather than level)
   */
  ScenarioGenerationBase(const MeanCov &norm,
			 arma::mat cone_A,
			 double alpha);

  /**
   * Constructor
   *
   * \param norm object representing Normal parameters of returns
   *             distribution
   * \param cone_A column major matrix of cone generators for feasible region.
   *               Each column represents a generator
   * \param num_gen number of cone generators
   * \param alpha quantile value (rather than level)
   */

  ScenarioGenerationBase(const MeanCov &norm,
			 double *cone_A, int num_gen,
			 double alpha);

  /**
   * Copy constructor
   */
  ScenarioGenerationBase(const ScenarioGenerationBase &);

  /**
   * Destructor
   */
  ~ScenarioGenerationBase();

  /**
   * Calculates for the transformation of a specified
   * point the distance from the transformed finite cone
   */
  double projected_distance(arma::vec point) const;

  /**
   * Tests whether a point is in the region of aggregation
   */
  bool in_agg_region(arma::vec point) const;

  /**
   * Prints a string representation to a specified output stream.
   */
  void pprint(std::ostream& = std::cout) const;

protected:
  int dim;
  int num_gen;
private:
  arma::vec mean;
  arma::mat P;
  arma::mat inv_P;
  ConeProjector *proj;
  double alpha;
};
    
#endif
