#include "AggregationSampling.hpp"
#include "julia_wrapper.hpp"

/**
 * Reduces scenario set by aggregating all scenarios in elliptical non-risk region.
 *
 * \param[out] new scenarios pointer to output scenarios (dim x num_scen preallocated double array)
 * \param[out] new_probs pointer to output probabilities (num_scen preallocated double array)
 * \param[in] old scenarios pointer to input scenarios
 * \param[in] new_probs pointer to output probabilities
 * \param dim dimension of problem
 * \param num_scen number of scenarios to be generated
 * \param num_gen number of generators for cone representing feasible region
 * \param mean mean vector used for elliptical risk region
 * \param cov covariance matrix used for elliptical risk region
 * \param cone_A matrix whose rows are the cone generators of feasible region
 * \param alpha cut-off distance for elliptical risk region
 */
int aggregation_julia(double *new_scenarios, double *new_probs, 
		      double *old_scenarios, double *old_probs, 
		      int dim, int num_scen, int num_gen,
		      double *mean, double *cov, double *cone_A,
		      double alpha)
{
  MeanCov meanCov(dim, mean, cov);
  AggregationSampling agg(meanCov, cone_A, num_gen, alpha);
  
  ScenarioSet scen(old_scenarios, old_probs, dim, num_scen);
  ScenarioSet new_scen = agg(scen);
  
  int new_num_scen  = new_scen.get_num_scen();
  for (int s = 0; s < new_num_scen; ++s) {
    new_probs[s] = new_scen.get_probability(s);
    for (int i = 0; i < dim; ++i)
      new_scenarios[s*dim + i] = new_scen.get_scenario(s,i);
  }
  return new_num_scen;
}


// /**
//  * Aggregation sampling using specified sampler.
//  * To use Julia callback for sampler, we use the "passthrough" idiom.
//  *
//  * \param[out] scenarios pointer to output scenarios (dim x num_scen preallocated double array)
//  * \param[out] probs pointer to output probabilities (num_scen preallocated double array)
//  * \param dim dimension of problem
//  * \param num_scen number of scenarios to be generated
//  * \param num_gen number of generators for cone representing feasible region
//  * \param mean mean vector used for elliptical risk region
//  * \param cov covariance matrix used for elliptical risk region
//  * \param cone_A matrix whose rows are the cone generators of feasible region
//  * \param alpha cut-off distance for elliptical risk region
//  * \param thunk parameter used for Julia passthrough function
//  * \param sampler sampling function (with pass-through parameter)
//  */
// int aggregation_sampling_julia(double *scenarios, double *probs,
// 			       int dim, int num_scen, int num_gen,
// 			       double *mean, double, *cov, double *cone_A,
// 			       double alpha, void *thunk,
// 			       void (*sampler)(void *thunk, double *scen))
// {
//   MeanCov meanCov(dim, mean, cov);
//   AggregationSampling agg(meanCov, cone_A, num_gen, alpha);
  
//   ScenarioSet scenarios = agg(num_scen);
  
  
// }
