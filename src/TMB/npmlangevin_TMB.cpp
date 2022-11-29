#define TMB_LIB_INIT R_init_npmlangevin
#include <TMB.hpp>
using namespace density;

#include "include/covariance.hpp"

#include "model/covariance_1d_deriv.hpp"
#include "other/covariance_exploration.hpp"


template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if( model == "covariance_exploration" ) {
    return covariance_exploration(this);
  } else if( model == "covariance_1d_deriv" ) {
    return covariance_1d_deriv(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}
