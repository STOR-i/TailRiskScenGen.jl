#include "MVNormal.hpp"
#include "AggregationSampling.hpp"

AggregationSampling::AggregationSampling(MeanCov &norm,
					       arma::mat cone_A,
					       double alpha)
  : ScenarioGenerationBase(norm, cone_A, alpha),
    sampler(new MVNormal(norm))
{}

AggregationSampling::AggregationSampling(MeanCov &norm,
					 double *cone_A, int num_gen,
					 double alpha)
  : ScenarioGenerationBase(norm, cone_A, num_gen, alpha),
    sampler(new MVNormal(norm))
{}


AggregationSampling::AggregationSampling(const AggregationSampling &other) 
  : ScenarioGenerationBase(other),
    sampler(new MVNormal(*(other.sampler)))
{}

AggregationSampling::~AggregationSampling()
{
  delete sampler;
}

ScenarioSet AggregationSampling::operator()(int num_scen)
{
  int non_agg_counter = 0;
  int agg_counter = 0;
  arma::mat scenarios(dim, num_scen);
  arma::vec probs(num_scen);
  arma::vec y(dim);
  arma::vec sum = arma::vec(dim, arma::fill::zeros);
  // Pointer stuff
  double *ptr = scenarios.memptr();
  
  while (non_agg_counter < num_scen - 1) {
    sampler->operator()(ptr);
    y = arma::vec(ptr, dim, false, true);

    if (in_agg_region(y)) {
      sum += y;
      ++agg_counter;
    }
    else {
      ++non_agg_counter;
      ptr += dim;
    }
  }
  if (agg_counter != 0) {
    scenarios.col(num_scen - 1) = sum/agg_counter;
    for (int s = 0; s < num_scen-1; ++s) 
      probs[s] = 1.0/(agg_counter + non_agg_counter);
    probs[num_scen-1] = float(agg_counter)/(agg_counter + non_agg_counter);
  }
  else {
    sampler->operator()(ptr);
    ++non_agg_counter;
    for (int s = 0; s < num_scen; ++s) probs[s] = 1.0/num_scen;
  }
  return ScenarioSet(scenarios, probs);
}

ScenarioSet AggregationSampling::operator()(const ScenarioSet &scen_set){
  int agg_counter = 0, non_agg_counter = 0;
  int orig_num_scen = scen_set.get_num_scen();
  int dim = scen_set.get_dim();
  arma::mat scenarios(dim, orig_num_scen);
  arma::vec probs(orig_num_scen);
  double agg_prob = 0;
  arma::vec y;
  arma::vec sum = arma::vec(dim, arma::fill::zeros);
  
  for (int s = 0; s < orig_num_scen; ++s) {
    y = scen_set.get_scenario(s);
    
    if (in_agg_region(y)) {
      sum += y;
      agg_prob += scen_set.get_probability(s);
      ++agg_counter;
    }
    else {
      scenarios.col(non_agg_counter) = y;
      probs[non_agg_counter] = scen_set.get_probability(s);
      ++non_agg_counter;
    }
  }
  if (agg_counter != 0) {
    scenarios.col(non_agg_counter) = sum/agg_counter;
    probs[non_agg_counter] = agg_prob;
    scenarios.resize(dim, non_agg_counter+1);
    probs.resize(non_agg_counter+1);
  }

  return ScenarioSet(scenarios, probs);
}
