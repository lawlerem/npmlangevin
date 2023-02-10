#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
matrix<Type> conjugate(matrix<Type> Sigma, matrix<Type> K, int q) {
  matrix<Type> diagK(2, 2);
  diagK << K(q, 0), 0.0, 0.0, K(q, 1);
  return diagK * Sigma * diagK;
}

template<class Type>
Type starve_npmlangevin(objective_function<Type>* obj) {
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
  REPORT(w);


  // Spatial field for utilization / gradient
  nngp<Type> field(g, w, boundary, cv);
  Type field_ll = field.loglikelihood();


  // Predictions for spatial field
  DATA_STRUCT(pwg, pred_graph);
  vector<Type> pw(pwg.var.size());

  for(int i = 0; i < pw.size(); i++) {
    pw(i) = field.predict(
      pwg.var(i),
      pwg.coord.row(i),
      pwg.parents(i)
    );
  }
  REPORT(pw);
  ADREPORT(pw);


  // Movement Process
  DATA_MATRIX(coordinates);
  DATA_STRUCT(field_neighbours, vmint);

  matrix<Type> coord_gradients(coordinates.rows(), coordinates.cols());

  for(int t = 0; t < coord_gradients.rows(); t++) {
    for(int v = 0; v < coord_gradients.cols(); v++) {
      coord_gradients(t, v) = field.predict(
        v + 1, // v = 0 = util distr, 1 = dx, 1 = dy
        vector<Type>(coordinates.row(t)),
        field_neighbours.x(t)
      );
    }
  }
  REPORT(coord_gradients);

  PARAMETER_MATRIX(random_walk);
  Type rw_ll = 0.0;
  for(int t = 0; t < random_walk.rows(); t++) {
    for(int v = 0; v < random_walk.cols(); v++) {
      rw_ll += dnorm(random_walk(t, v), Type(0.0), Type(1.0), true);
    }
  }

  DATA_VECTOR(time);
  PARAMETER(log_gamma);
  Type gamma = exp(log_gamma);
  ADREPORT(gamma);

  matrix<Type> location_difference_means(random_walk.rows(), random_walk.cols());
  for(int t = 0; t < location_difference_means.rows(); t++) {
    for(int v = 0; v < location_difference_means.cols(); v++) {
      location_difference_means(t, v) = 0.5 * (time(t + 1) - time(t)) * coord_gradients(t, v) + gamma * sqrt(time(t + 1) - time(t)) * random_walk(t, v);
    }
  }

  // Observation Noise + random walk movement noise
  DATA_MATRIX(location_differences);
  DATA_IVECTOR(location_quality_class);

  DATA_MATRIX(K);
  PARAMETER_VECTOR(working_ping_cov_pars);
  vector<Type> ping_cov_pars = exp(working_ping_cov_pars);
  ping_cov_pars(1) = 2 * invlogit(working_ping_cov_pars(1)) - 1.0;

  matrix<Type> ping_cov(2, 2);
  ping_cov << pow(ping_cov_pars(0), 2), ping_cov_pars(1) * ping_cov_pars(0) * ping_cov_pars(2),
    ping_cov_pars(1) * ping_cov_pars(0) * ping_cov_pars(2), pow(ping_cov_pars(2), 2);
  ADREPORT(ping_cov);

  Type ping_ll = 0.0;
  for(int t = 0; t < location_differences.rows(); t++) {
    matrix<Type> Sigma = conjugate(ping_cov, K, location_quality_class(t)) + conjugate(ping_cov, K, location_quality_class(t + 1));
    ping_ll -= MVNORM<Type>(Sigma)(vector<Type>(location_differences.row(t)));
  }

return -1.0 * (field_ll + rw_ll + ping_ll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this