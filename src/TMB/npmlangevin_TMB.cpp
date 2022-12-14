#define TMB_LIB_INIT R_init_npmlangevin
#include <TMB.hpp>
using namespace density;

#include "include/utilities.hpp"
#include "include/boundary_mean.hpp"
#include "include/covariance.hpp"
#include "include/conditional_normal.hpp"
#include "include/graph.hpp"
#include "include/pred_graph.hpp"
#include "include/nngp.hpp"
#include "include/loc_track.hpp"
#include "include/loc_observations.hpp"


#include "model/covariance_1d_deriv.hpp"
#include "model/nngp_model.hpp"
#include "model/random_walk.hpp"
#include "model/langevin_diffusion.hpp"

#include "other/covariance_exploration.hpp"


template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if( model == "covariance_exploration" ) {
    return covariance_exploration(this);
  } else if( model == "covariance_1d_deriv" ) {
    return covariance_1d_deriv(this);
  } else if( model == "nngp_model" ) {
    return nngp_model(this);
  } else if( model == "random_walk" ) {
    return random_walk(this);
  } else if( model == "langevin_diffusion" ) {
    return langevin_diffusion(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}
