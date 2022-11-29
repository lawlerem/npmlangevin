#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type covariance_1d_deriv(objective_function<Type>* obj) {
  DATA_INTEGER(cv_code);
  DATA_VECTOR(cv_pars);

  PARAMETER(x);

  covariance<Type> f {cv_pars, cv_code};

  vector<Type> x1(1);
  x1 << 0.0;

  vector<Type> x2(1);
  x2 << x;

  return f.gradient(x1, x2)(1);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
