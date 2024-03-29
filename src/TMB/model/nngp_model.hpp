#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type nngp_model(objective_function<Type>* obj) {
  DATA_INTEGER(cv_code);
  DATA_STRUCT(g, nngp_graph);
  DATA_ARRAY(y);
  DATA_STRUCT(pwg, pred_graph);

  PARAMETER_VECTOR(boundary_x);
  PARAMETER_VECTOR(boundary_y);
  PARAMETER(working_boundary_sharpness);

  Type boundary_sharpness = exp(working_boundary_sharpness);
  ADREPORT(boundary_sharpness);

  boundary_mean<Type> boundary {boundary_x, boundary_y, boundary_sharpness};

  PARAMETER_VECTOR(working_cv_pars);
  PARAMETER_ARRAY(w);
  vector<Type> pw(pwg.var.size());

  vector<Type> cv_pars = exp(working_cv_pars);
  ADREPORT(cv_pars);

  covariance<Type> cv {cv_pars, cv_code};
  nngp<Type> field(g, w, boundary, cv);

  Type field_ll = field.loglikelihood();
  SIMULATE{
    w = field.simulate();
    REPORT(w);
  }

  Type obs_ll = 0.0;
  for(int i = 0; i < y.dim(0); i++) {
    for(int j = 0; j < y.dim(1); j++) {
      for(int k = 0; k < y.dim(2); k++) {
        obs_ll += dnorm(y(i, j, k), w(i, j, k), Type(1.0), true);
      }
    }
  }
  SIMULATE{
    for(int i = 0; i < y.dim(0); i++) {
      for(int j = 0; j < y.dim(1); j++) {
        for(int k = 0; k < y.dim(2); k++) {
          y(i, j, k) = rnorm(w(i, j, k), Type(1.0));
        }
      }
    }
    REPORT(y);
  }

  for(int i = 0; i < pw.size(); i++) {
    pw(i) = field.predict(
      pwg.var(i),
      pwg.coord.row(i),
      pwg.parents(i)
    );
  }
  ADREPORT(pw);

  Type ll = field_ll + obs_ll;
  REPORT(ll);
  REPORT(field_ll);
  REPORT(obs_ll);

  return -1.0 * ll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
