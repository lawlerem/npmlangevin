#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type langevin_diffusion(objective_function<Type>* obj) {
  // Boundary effects
  PARAMETER_VECTOR(boundary_x);
  PARAMETER_VECTOR(boundary_y);
  PARAMETER(working_boundary_sharpness);
  Type boundary_sharpness = exp(working_boundary_sharpness);
  ADREPORT(boundary_sharpness);
  boundary_mean<Type> boundary {boundary_x, boundary_y, boundary_sharpness};

  // Covariance function
  DATA_INTEGER(cv_code);
  PARAMETER_VECTOR(working_cv_pars);
  vector<Type> cv_pars = exp(working_cv_pars);
  ADREPORT(cv_pars);
  covariance<Type> cv {cv_pars, cv_code};

  // Nearest neighbour graph and random effects
  DATA_STRUCT(g, nngp_graph);
  PARAMETER_ARRAY(w);

  // Spatial field for utilization / gradient
  nngp<Type> field(g, w, boundary, cv);

  Type field_ll = field.loglikelihood();
  SIMULATE{
    w = field.simulate();
    REPORT(w);
  }

  // Predictions for spatial field
  DATA_STRUCT(pwg, pred_graph);
  PARAMETER_VECTOR(pw);
  Type pred_ll = 0.0;
  for(int i = 0; i < pw.size(); i++) {
    field.predict(
      pw(i),
      pwg.var(i),
      pwg.coord.row(i),
      pwg.parents(i),
      pred_ll
    );
  }

  SIMULATE{
    for(int i = 0; i < pw.size(); i++) {
      pw(i) = field.simulate_predict(
        pwg.var(i),
        pwg.coord.row(i),
        pwg.parents(i),
        true
      );
    }
    REPORT(pw);
  }

  // True Movement Path
  PARAMETER_MATRIX(true_coord);
  DATA_STRUCT(field_neighbours, vmint);
  DATA_VECTOR(true_time);
  PARAMETER(log_gamma);
  Type gamma = exp(log_gamma);
  ADREPORT(gamma);

  loc_track<Type> track {true_coord, field_neighbours, true_time, gamma};
  Type track_ll = track.loglikelihood(field);

  SIMULATE{
    true_coord = track.simulate(field);
    REPORT(true_coord);
  }

  // Observed Locations
  DATA_STRUCT(pings, loc_observations);

  PARAMETER_VECTOR(working_ping_cov_pars);
  vector<Type> ping_cov_pars = exp(working_ping_cov_pars);
  ping_cov_pars(1) = 2 * invlogit(working_ping_cov_pars(1)) - 1.0;
  matrix<Type> ping_cov(2, 2);
  ping_cov << pow(ping_cov_pars(0), 2), ping_cov_pars(1) * ping_cov_pars(0) * ping_cov_pars(2),
    ping_cov_pars(1) * ping_cov_pars(0) * ping_cov_pars(2), pow(ping_cov_pars(2), 2);
  ADREPORT(ping_cov);

  Type pings_ll = pings.loglikelihood(ping_cov, track);
  SIMULATE{
    matrix<Type> sim_pings = pings.simulate(ping_cov, track);
    REPORT(sim_pings);
  }

  REPORT(field_ll);
  REPORT(pred_ll);
  REPORT(track_ll);
  REPORT(pings_ll);

  return -1.0 * (field_ll + pred_ll + track_ll + pings_ll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this