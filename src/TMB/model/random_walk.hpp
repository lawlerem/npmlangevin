#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type random_walk(objective_function<Type>* obj) {
  PARAMETER_MATRIX(true_loc);
  DATA_VECTOR(true_time);
  PARAMETER(log_gamma);

  Type gamma = exp(log_gamma);
  ADREPORT(gamma);

  loc_track<Type> true_track {true_loc, true_time, gamma};


  DATA_STRUCT(loc_obs, loc_observations);

  PARAMETER_VECTOR(working_obs_cov_pars);
  vector<Type> obs_cov_pars = exp(working_obs_cov_pars);
  obs_cov_pars(1) = 2 * invlogit(working_obs_cov_pars(1)) - 1.0;
  matrix<Type> Sigma(2, 2);
  Sigma << pow(obs_cov_pars(0), 2), obs_cov_pars(1) * obs_cov_pars(0) * obs_cov_pars(2),
    obs_cov_pars(1) * obs_cov_pars(0) * obs_cov_pars(2), pow(obs_cov_pars(2), 2);

  ADREPORT(Sigma);

  Type proc_ll = true_track.loglikelihood();
  Type obs_ll = loc_obs.loglikelihood(Sigma, true_track);
  Type ll = proc_ll + obs_ll;

  REPORT(ll);
  REPORT(proc_ll);
  REPORT(obs_ll);

  SIMULATE{
    true_loc = true_track.simulate();
    matrix<Type> sim_obs = loc_obs.simulate(Sigma, true_track);
    REPORT(true_loc);
    REPORT(sim_obs);
  }

  return -1.0 * ll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
