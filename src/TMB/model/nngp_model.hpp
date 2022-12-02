#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type nngp_model(objective_function<Type>* obj) {
  DATA_INTEGER(cv_code);
  DATA_STRUCT(g, nngp_graph);
  DATA_ARRAY(y);

  PARAMETER_VECTOR(working_cv_pars);
  PARAMETER_ARRAY(w);

  vector<Type> cv_pars = exp(working_cv_pars);
  ADREPORT(cv_pars);

  covariance<Type> cv {cv_pars, cv_code};
  nngp<Type> field(g, w, cv);

  Type nll = -1.0 * field.loglikelihood();
  SIMULATE{
    w = field.simulate();
    REPORT(w);
  }

  for(int i = 0; i < y.dim(0); i++) {
    for(int j = 0; j < y.dim(1); j++) {
      for(int k = 0; k < y.dim(2); k++) {
        nll -= dnorm(y(i, j, k), w(i, j, k), Type(1.0), true);
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

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
