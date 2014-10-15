#include "AggregationSampling.hpp"
#include "julia_wrapper.hpp"

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
