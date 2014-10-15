#ifndef AGGREGATION_SAMPLING_HPP
#define AGGREGATION_SAMPLING_HPP

#include "AggregationSampling.hpp"

int aggregation_julia(double *new_scenarios, double *old_scenarios, int dim,
			int num_scen, int num_gen,
			double *mean, double *cov, double *cone_A,
			double alpha);

#endif
