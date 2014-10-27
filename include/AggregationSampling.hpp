#ifndef AGGREGATION_SAMPLING_HPP
#define AGGREGATION_SAMPLING_HPP

#include <functional>
#include "MeanCov.hpp"
#include "ScenarioGenerationBase.hpp"
#include "ScenarioSet.hpp"

/**
 * Provides methods for simple scneario aggregation
 * where the aggregation region is represented with one point.
 */
class AggregationSampling : public ScenarioGenerationBase {
public:
  /**
   * Constructor
   */
  AggregationSampling(MeanCov &norm,
		      arma::mat cone_A,
		      double alpha);

  AggregationSampling(MeanCov &norm,
		      double *cone_A, int num_gen,
		      double alpha);

  AggregationSampling(std::function< void(double*) > sampler,
		      MeanCov &norm,
		      double *cone_A, int num_gen,
		      double alpha);

  /**
   * Copy constructor
   */
  AggregationSampling(const AggregationSampling &);

  /**
   * Destructor
   */
  ~AggregationSampling();

  /**
   * Constructs a scenario set by sampling from the specified
   * Normal distribution, aggregating points in aggregation
   * region and storing other points, until we have the required
   * number of scenarios.
   *
   * \param num_scen number of scenarios in constructed set
   * \return generated scenario set
   */
  ScenarioSet operator()(int num_scen);

  /**
   * Reduces a given scenario set by aggregating
   * any scenarios in the region of aggregation
   *
   * \param scen_set original scenario set
   * \return reduced scenario set
   */
  ScenarioSet operator()(const ScenarioSet &scen_set);
private:
  std::function< void(double *) > sampler;
};

#endif
