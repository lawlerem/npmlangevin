#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type nngp_model(objective_function<Type>* obj) {
  DATA_INTEGER(cv_code);
  DATA_STRUCT(g, nngp_graph);
  DATA_ARRAY(y);
  DATA_STRUCT(pwg, pred_graph);

  PARAMETER_VECTOR(working_cv_pars);
  PARAMETER_ARRAY(w);
  PARAMETER_VECTOR(pw);

  vector<Type> cv_pars = exp(working_cv_pars);
  ADREPORT(cv_pars);

  covariance<Type> cv {cv_pars, cv_code};
  nngp<Type> field(g, w, cv);

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

  Type ll = field_ll + obs_ll + pred_ll;
  REPORT(ll);
  REPORT(field_ll);
  REPORT(obs_ll);
  REPORT(pred_ll);

  return -1.0 * ll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
