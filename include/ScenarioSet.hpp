#ifndef SCENARIO_SET_HPP
#define SCENARIO_SET_HPP

#include <armadillo>
#include "MeanCov.hpp"


/**
 * Class for representing general scenario sets.
 */

class ScenarioSet {
public:
  /**
   * Default constructor.
   * Construct empty scenario set
   */
  ScenarioSet() {}

  /**
   * Constructs scenario set from a matrix of scenarios
   * and a vector of associated probabilities.
   * 
   * Input is checked for consistency.
   *
   * \param scenarios matrix of scenarios where
   *         each column represents one scenario
   * \param vector of probabilities
   */
  ScenarioSet(arma::mat scenarios, arma::vec probs);

  /**
   * Constructs scenario set from C pointer arrays
   *
   * \param scenarios column-major C-array containing scenarios
   *                  Each column is a scenario
   */
  ScenarioSet(double *scenarios, double *probs, int dim, int num_scen);
 
  /**
   * Constructs scenario set with equal probabilities.
   * \param scenarios matrix of scenarios where
   *        each column represents one scenario
   */
  ScenarioSet(arma::mat scenarios);
  
  int get_dim() const {return dim;};

  int get_num_scen() const {return num_scen;};

  /**
   * Get the value of an ordinate of a particular scenario.
   * \param s index of scenario
   * \param i index of ordinate
   */
  double get_scenario(int s, int i) const {return scenarios(i,s);}

  /**
   * Get a specified scenario.
   * \param s index of scenario
   */
  arma::vec get_scenario(int s) const {return scenarios.col(s);}

  
  /**
   * Get probability of specified scenario.
   * \param s index of scenario
   */
  double get_probability(int s) const {return probs(s);}

  /**
   * Print scenario set to a specified output stream
   */
  void pprint(std::ostream &out = std::cout) const;

  /**
   * Clear current scenario set and read new set from specified file
   * \param fn file name
   */
  void read(std::string fn);

private:
  arma::mat scenarios;
  arma::vec probs;
  int dim;
  int num_scen;
};

/*
 * Function which constructs a random
 * scenario set by sampling normal variates
 */
ScenarioSet rnorm_ScenarioSet(const MeanCov &norm, int num_scen);

#endif
